# RNA-seq Analysis Web Service

A web-based RNA-seq analysis pipeline powered by Nextflow, FastAPI, and React.

## Quick Start

### Prerequisites

- Docker and Docker Compose installed
- At least 8GB RAM available for Docker
- Reference genome FASTA file (*.genome.fa) in the workspace root
- GTF annotation file (optional, or will use default)

### Starting the Service

1. Clone the repository and navigate to the workspace directory

2. Copy the example environment file:
```bash
cp .env.example .env
```

3. Edit `.env` and change the JWT_SECRET to a secure random string

4. Start all services using Docker Compose:
```bash
docker-compose up -d
```

5. Access the web interface at http://localhost:3000

6. The API documentation is available at http://localhost:8000/docs

### Stopping the Service

```bash
docker-compose down
```

To also remove data volumes:
```bash
docker-compose down -v
```

## Using the Web Interface

### First Time Setup

1. Open http://localhost:3000 in your browser
2. Click on the "Register" tab
3. Create a new account with username, email, and password
4. You'll be automatically logged in after registration

### Submitting a Job

1. Click "New Job" in the navigation bar
2. Select the read mode:
   - **Paired-end**: For paired-end FASTQ files (default)
   - **Single-end**: For single-end FASTQ files
3. (Optional) Specify custom GTF file path or STAR index path
4. Drag and drop your FASTQ files or click to select them
5. Click "Create and Start Job"

### Monitoring Jobs

1. The dashboard shows all your jobs with their current status
2. Click the eye icon to view detailed job information
3. The job monitor shows:
   - Real-time status updates via WebSocket
   - Live Nextflow log output
   - Results files when the job completes

### Downloading Results

1. Navigate to a completed job
2. Click the "Results" tab
3. Click "Download" next to any result file

## Architecture

The service consists of three main components:

### Backend (FastAPI)
- REST API for job management
- WebSocket for real-time updates
- SQLite database for user and job data
- Executes Nextflow pipeline as subprocess
- Located in `./backend/`

### Frontend (React)
- Material-UI based responsive interface
- File upload with drag-and-drop
- Real-time job monitoring
- Authentication and authorization
- Located in `./frontend/`

### Pipeline (Nextflow)
- Existing RNA-seq analysis workflow
- Runs in Docker containers
- Stages: TRIM → STAR_INDEX → ALIGN → FEATURECOUNT → MULTIQC
- Uses existing modules in `./modules/`

## API Endpoints

### Authentication
- `POST /api/auth/register` - Create new user account
- `POST /api/auth/login` - Login and get JWT token

### Jobs
- `POST /api/jobs` - Create a new job
- `POST /api/jobs/{id}/upload` - Upload FASTQ file to job
- `POST /api/jobs/{id}/start` - Start job execution
- `GET /api/jobs` - List all user's jobs
- `GET /api/jobs/{id}` - Get job details
- `GET /api/jobs/{id}/log` - Get job log output
- `GET /api/jobs/{id}/results` - List result files
- `GET /api/jobs/{id}/results/{path}` - Download result file
- `DELETE /api/jobs/{id}` - Cancel running job

### WebSocket
- `WS /ws/jobs/{id}` - Real-time job status updates

## Configuration

### Environment Variables

Create a `.env` file in the workspace root:

```bash
JWT_SECRET=your-secure-secret-key-here
DATABASE_URL=sqlite:///./rnaseq.db
REACT_APP_API_URL=http://localhost:8000
```

### Pipeline Parameters

Default parameters can be modified in `nextflow.config`:
- CPU and memory allocation
- Docker containers for each stage
- Default paths for reference files

### Job Storage

- Job files are stored in `/workspace/jobs/{job_id}/`
- Each job has its own directory structure:
  - `data/` - Uploaded FASTQ files
  - `results/` - Pipeline output files
  - `work/` - Nextflow working directory (temporary)

## Development

### Backend Development

```bash
cd backend
pip install -r requirements.txt
uvicorn app.main:app --reload
```

### Frontend Development

```bash
cd frontend
npm install
npm start
```

### Running Tests

```bash
cd tests
bash run_tests.sh
```

## Troubleshooting

### Backend won't start
- Check if port 8000 is available
- Ensure Docker socket is accessible
- Check logs: `docker-compose logs backend`

### Frontend won't start
- Check if port 3000 is available
- Verify backend is running
- Check logs: `docker-compose logs frontend`

### Pipeline fails
- Ensure reference genome file exists
- Check available memory (STAR requires 30-40GB)
- Review Nextflow logs in job monitor
- Verify Docker containers can be pulled

### WebSocket connection fails
- Check if both backend and frontend are running
- Verify no firewall blocking WebSocket connections
- Check browser console for errors

## Security Notes

For production deployment:
- Change JWT_SECRET to a strong random string
- Use PostgreSQL instead of SQLite
- Enable HTTPS/TLS
- Implement rate limiting
- Add proper authentication middleware
- Configure CORS properly
- Use environment-specific .env files

## Support

For issues related to:
- Web service: Check Docker logs and browser console
- Pipeline: See original README.md and Nextflow documentation
- Nextflow configuration: See nextflow.config comments
