from fastapi import FastAPI, Depends, HTTPException, status, UploadFile, File, WebSocket, WebSocketDisconnect, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from typing import List, Optional
from pydantic import BaseModel, EmailStr, validator
from datetime import timedelta
import os
import asyncio
from pathlib import Path

from . import models, auth
from .database import get_db, init_db
from .pipeline import PipelineExecutor, JOBS_DIR
from .middleware.security import SecurityHeadersMiddleware
from .middleware.rate_limit import limiter, get_limiter
from .logging_config import logger
from slowapi.errors import RateLimitExceeded
from slowapi import _rate_limit_exceeded_handler
from prometheus_fastapi_instrumentator import Instrumentator

app = FastAPI(title="RNA-seq Analysis API", version="1.0.0")

Instrumentator().instrument(app).expose(app, endpoint="/metrics")

app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

allowed_origins = os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
)

app.add_middleware(SecurityHeadersMiddleware)

@app.on_event("startup")
async def startup_event():
    logger.info("Starting RNA-seq Analysis API")
    init_db()
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("Application startup complete")

class UserCreate(BaseModel):
    username: str
    email: EmailStr
    password: str

class UserLogin(BaseModel):
    username: str
    password: str

class Token(BaseModel):
    access_token: str
    token_type: str

class JobCreate(BaseModel):
    single_end: bool = False
    gtf: Optional[str] = None
    star_index: Optional[str] = None
    
    @validator('gtf', 'star_index')
    def validate_path(cls, v):
        if v is not None:
            if '..' in v or v.startswith('/'):
                raise ValueError('Invalid path: directory traversal not allowed')
        return v

class JobResponse(BaseModel):
    id: int
    status: str
    parameters: dict
    created_at: str
    updated_at: str
    results_path: Optional[str]
    error_message: Optional[str]
    
    class Config:
        from_attributes = True

@app.post("/api/auth/register", response_model=Token)
@limiter.limit("5/minute")
def register(request: Request, user: UserCreate, db: Session = Depends(get_db)):
    logger.info(f"Registration attempt for username: {user.username}")
    
    existing_user = db.query(models.User).filter(
        (models.User.username == user.username) | (models.User.email == user.email)
    ).first()
    
    if existing_user:
        logger.warning(f"Registration failed: username/email already exists - {user.username}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Username or email already registered"
        )
    
    hashed_password = auth.get_password_hash(user.password)
    db_user = models.User(
        username=user.username,
        email=user.email,
        hashed_password=hashed_password
    )
    db.add(db_user)
    db.commit()
    db.refresh(db_user)
    
    logger.info(f"User registered successfully: {user.username}")
    
    access_token = auth.create_access_token(data={"sub": user.username})
    return {"access_token": access_token, "token_type": "bearer"}

@app.post("/api/auth/login", response_model=Token)
@limiter.limit("5/minute")
def login(request: Request, user: UserLogin, db: Session = Depends(get_db)):
    db_user = db.query(models.User).filter(models.User.username == user.username).first()
    
    if not db_user or not auth.verify_password(user.password, db_user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password"
        )
    
    access_token = auth.create_access_token(data={"sub": user.username})
    return {"access_token": access_token, "token_type": "bearer"}

@app.post("/api/jobs", response_model=JobResponse)
@limiter.limit("10/hour")
async def create_job(
    request: Request,
    job_data: JobCreate,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = models.Job(
        user_id=current_user.id,
        status="pending",
        parameters=job_data.dict()
    )
    db.add(job)
    db.commit()
    db.refresh(job)
    
    return JobResponse(
        id=job.id,
        status=job.status,
        parameters=job.parameters,
        created_at=job.created_at.isoformat(),
        updated_at=job.updated_at.isoformat(),
        results_path=job.results_path,
        error_message=job.error_message
    )

@app.post("/api/jobs/{job_id}/upload")
@limiter.limit("20/minute")
async def upload_file(
    request: Request,
    job_id: int,
    file: UploadFile = File(...),
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    MAX_FILE_SIZE = 10 * 1024 * 1024 * 1024
    
    if not file.filename.endswith(('.fastq.gz', '.fq.gz')):
        raise HTTPException(
            status_code=400,
            detail="Only FASTQ files (.fastq.gz or .fq.gz) are allowed"
        )
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job_dir = JOBS_DIR / str(job_id) / "data"
    job_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = job_dir / file.filename
    with open(file_path, "wb") as f:
        content = await file.read()
        f.write(content)
    
    return {"filename": file.filename, "size": len(content)}

@app.post("/api/jobs/{job_id}/start")
async def start_job(
    job_id: int,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if job.status != "pending":
        raise HTTPException(status_code=400, detail="Job already started")
    
    from .tasks import run_pipeline_task
    task = run_pipeline_task.delay(job_id)
    
    return {"message": "Job queued", "job_id": job_id, "task_id": task.id}

@app.get("/api/jobs", response_model=List[JobResponse])
def list_jobs(
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    jobs = db.query(models.Job).filter(
        models.Job.user_id == current_user.id
    ).order_by(models.Job.created_at.desc()).all()
    
    return [
        JobResponse(
            id=job.id,
            status=job.status,
            parameters=job.parameters,
            created_at=job.created_at.isoformat(),
            updated_at=job.updated_at.isoformat(),
            results_path=job.results_path,
            error_message=job.error_message
        )
        for job in jobs
    ]

@app.get("/api/jobs/{job_id}", response_model=JobResponse)
def get_job(
    job_id: int,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return JobResponse(
        id=job.id,
        status=job.status,
        parameters=job.parameters,
        created_at=job.created_at.isoformat(),
        updated_at=job.updated_at.isoformat(),
        results_path=job.results_path,
        error_message=job.error_message
    )

@app.get("/api/jobs/{job_id}/log")
def get_job_log(
    job_id: int,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return {"log": job.nextflow_log or ""}

@app.get("/api/jobs/{job_id}/results")
def list_results(
    job_id: int,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job or not job.results_path:
        raise HTTPException(status_code=404, detail="Results not found")
    
    results_dir = Path(job.results_path)
    if not results_dir.exists():
        raise HTTPException(status_code=404, detail="Results directory not found")
    
    files = []
    for item in results_dir.rglob("*"):
        if item.is_file():
            rel_path = item.relative_to(results_dir)
            files.append({
                "path": str(rel_path),
                "size": item.stat().st_size,
                "modified": item.stat().st_mtime
            })
    
    return {"files": files}

@app.get("/api/jobs/{job_id}/results/{file_path:path}")
def download_result(
    job_id: int,
    file_path: str,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job or not job.results_path:
        raise HTTPException(status_code=404, detail="Results not found")
    
    file_full_path = Path(job.results_path) / file_path
    if not file_full_path.exists() or not file_full_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    
    return FileResponse(file_full_path)

@app.delete("/api/jobs/{job_id}")
def cancel_job(
    job_id: int,
    current_user: models.User = Depends(auth.get_current_user),
    db: Session = Depends(get_db)
):
    job = db.query(models.Job).filter(
        models.Job.id == job_id,
        models.Job.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if job.status == "running":
        executor = PipelineExecutor(job_id, db)
        executor.cancel_job()
    
    return {"message": "Job cancelled"}

@app.websocket("/ws/jobs/{job_id}")
async def websocket_job_status(websocket: WebSocket, job_id: int, db: Session = Depends(get_db)):
    await websocket.accept()
    
    try:
        while True:
            job = db.query(models.Job).filter(models.Job.id == job_id).first()
            if job:
                await websocket.send_json({
                    "status": job.status,
                    "updated_at": job.updated_at.isoformat(),
                    "log": job.nextflow_log or ""
                })
            
            if job and job.status in ["completed", "failed", "cancelled"]:
                break
            
            await asyncio.sleep(2)
            
    except WebSocketDisconnect:
        pass

@app.get("/")
def root():
    return {"message": "RNA-seq Analysis API", "version": "1.0.0"}

@app.get("/health")
def health(db: Session = Depends(get_db)):
    import shutil
    from sqlalchemy import text
    
    health_status = {
        "status": "healthy",
        "checks": {}
    }
    
    try:
        db.execute(text("SELECT 1"))
        health_status["checks"]["database"] = "healthy"
    except Exception as e:
        health_status["status"] = "degraded"
        health_status["checks"]["database"] = f"unhealthy: {str(e)}"
        logger.error(f"Database health check failed: {e}")
    
    try:
        disk_usage = shutil.disk_usage("/workspace/jobs")
        free_gb = disk_usage.free / (1024**3)
        if free_gb < 10:
            health_status["status"] = "degraded"
            health_status["checks"]["disk_space"] = f"low: {free_gb:.2f}GB free"
        else:
            health_status["checks"]["disk_space"] = f"healthy: {free_gb:.2f}GB free"
    except Exception as e:
        health_status["checks"]["disk_space"] = f"unknown: {str(e)}"
    
    try:
        import redis
        redis_client = redis.from_url(os.getenv("REDIS_URL", "redis://localhost:6379/0"))
        redis_client.ping()
        health_status["checks"]["redis"] = "healthy"
    except Exception as e:
        health_status["status"] = "degraded"
        health_status["checks"]["redis"] = f"unhealthy: {str(e)}"
        logger.error(f"Redis health check failed: {e}")
    
    return health_status
