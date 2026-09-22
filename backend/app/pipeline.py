import asyncio
import os
import shutil
from pathlib import Path
from typing import Optional
from sqlalchemy.orm import Session
from . import models

JOBS_DIR = Path("/workspace/jobs")
NEXTFLOW_BIN = "nextflow"
PIPELINE_DIR = Path("/workspace")

class PipelineExecutor:
    def __init__(self, job_id: int, db: Session):
        self.job_id = job_id
        self.db = db
        self.job_dir = JOBS_DIR / str(job_id)
        self.data_dir = self.job_dir / "data"
        self.results_dir = self.job_dir / "results"
        self.work_dir = self.job_dir / "work"
        
    def setup_directories(self):
        self.job_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)
        self.work_dir.mkdir(exist_ok=True)
        
    async def run_pipeline(self, parameters: dict):
        self.setup_directories()
        
        job = self.db.query(models.Job).filter(models.Job.id == self.job_id).first()
        if not job:
            raise ValueError(f"Job {self.job_id} not found")
        
        job.status = "running"
        self.db.commit()
        
        cmd = self._build_nextflow_command(parameters)
        
        try:
            process = await asyncio.create_subprocess_shell(
                cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
                cwd=str(PIPELINE_DIR)
            )
            
            log_lines = []
            async for line in process.stdout:
                line_str = line.decode('utf-8')
                log_lines.append(line_str)
                
                job.nextflow_log = ''.join(log_lines[-1000:])
                self.db.commit()
            
            await process.wait()
            
            if process.returncode == 0:
                job.status = "completed"
                job.results_path = str(self.results_dir)
            else:
                job.status = "failed"
                job.error_message = f"Pipeline failed with exit code {process.returncode}"
                
        except Exception as e:
            job.status = "failed"
            job.error_message = str(e)
        
        finally:
            self.db.commit()
            
    def _build_nextflow_command(self, parameters: dict) -> str:
        single_end = parameters.get('single_end', False)
        gtf = parameters.get('gtf')
        star_index = parameters.get('star_index')
        
        cmd_parts = [
            NEXTFLOW_BIN,
            "run",
            str(PIPELINE_DIR / "main.nf"),
            f"--outdir {self.results_dir}",
            f"--data_dir {self.data_dir}",
        ]
        
        if single_end:
            cmd_parts.append("--single_end true")
        
        if gtf:
            cmd_parts.append(f"--gtf {gtf}")
            
        if star_index:
            cmd_parts.append(f"--star_index {star_index}")
        
        cmd_parts.extend([
            f"-work-dir {self.work_dir}",
            "-with-docker"
        ])
        
        return " ".join(cmd_parts)
    
    def cancel_job(self):
        job = self.db.query(models.Job).filter(models.Job.id == self.job_id).first()
        if job and job.status == "running":
            job.status = "cancelled"
            self.db.commit()
            
    def cleanup(self):
        if self.work_dir.exists():
            shutil.rmtree(self.work_dir, ignore_errors=True)
