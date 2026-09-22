import asyncio
from celery import Task
from sqlalchemy.orm import Session
from .celery_app import celery_app
from .database import SessionLocal
from .pipeline import PipelineExecutor
from . import models


class DatabaseTask(Task):
    _db = None

    @property
    def db(self):
        if self._db is None:
            self._db = SessionLocal()
        return self._db

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()


@celery_app.task(bind=True, base=DatabaseTask, name="app.tasks.run_pipeline_task")
def run_pipeline_task(self, job_id: int):
    db = self.db
    
    job = db.query(models.Job).filter(models.Job.id == job_id).first()
    if not job:
        return {"error": f"Job {job_id} not found"}
    
    job.status = "running"
    db.commit()
    
    try:
        executor = PipelineExecutor(job_id, db)
        
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        loop.run_until_complete(executor.run_pipeline(job.parameters))
        loop.close()
        
        return {"job_id": job_id, "status": "completed"}
    
    except Exception as e:
        job.status = "failed"
        job.error_message = str(e)
        db.commit()
        
        return {"job_id": job_id, "status": "failed", "error": str(e)}


@celery_app.task(name="app.tasks.send_email_task")
def send_email_task(to_email: str, subject: str, body: str):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    import os
    
    smtp_host = os.getenv("SMTP_HOST")
    smtp_port = int(os.getenv("SMTP_PORT", "587"))
    smtp_user = os.getenv("SMTP_USER")
    smtp_password = os.getenv("SMTP_PASSWORD")
    
    if not all([smtp_host, smtp_user, smtp_password]):
        return {"status": "skipped", "reason": "SMTP not configured"}
    
    try:
        msg = MIMEMultipart()
        msg["From"] = smtp_user
        msg["To"] = to_email
        msg["Subject"] = subject
        msg.attach(MIMEText(body, "html"))
        
        server = smtplib.SMTP(smtp_host, smtp_port)
        server.starttls()
        server.login(smtp_user, smtp_password)
        server.send_message(msg)
        server.quit()
        
        return {"status": "sent", "to": to_email}
    
    except Exception as e:
        return {"status": "failed", "error": str(e)}


@celery_app.task(name="app.tasks.cleanup_old_jobs_task")
def cleanup_old_jobs_task(days: int = 30):
    from datetime import datetime, timedelta
    import shutil
    from pathlib import Path
    
    db = SessionLocal()
    try:
        cutoff_date = datetime.utcnow() - timedelta(days=days)
        
        old_jobs = db.query(models.Job).filter(
            models.Job.created_at < cutoff_date,
            models.Job.status.in_(["completed", "failed", "cancelled"])
        ).all()
        
        cleaned_count = 0
        for job in old_jobs:
            job_dir = Path(f"/workspace/jobs/{job.id}")
            if job_dir.exists():
                shutil.rmtree(job_dir, ignore_errors=True)
            
            db.delete(job)
            cleaned_count += 1
        
        db.commit()
        
        return {"cleaned": cleaned_count, "cutoff_date": cutoff_date.isoformat()}
    
    except Exception as e:
        db.rollback()
        return {"error": str(e)}
    
    finally:
        db.close()
