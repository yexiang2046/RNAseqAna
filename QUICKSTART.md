# Quick Start Guide

Get the RNA-seq web service running in 5 minutes!

## Prerequisites

- Docker and Docker Compose installed
- At least 8GB RAM available for Docker
- Reference genome file (*.genome.fa) in workspace root

## Step 1: Setup Environment

```bash
# Copy environment template
cp .env.example .env

# Edit and set a secure JWT secret
nano .env
# Change: JWT_SECRET=your-secure-random-string-here
```

## Step 2: Start Services

```bash
# Build and start all services
docker-compose up -d

# Check status
docker-compose ps

# View logs (optional)
docker-compose logs -f
```

Expected output:
```
NAME                  STATUS              PORTS
rnaseq-backend        Up (healthy)        0.0.0.0:8000->8000/tcp
rnaseq-frontend       Up (healthy)        0.0.0.0:3000->80/tcp
```

## Step 3: Access Web Interface

Open your browser to: **http://localhost:3000**

## Step 4: Create Account

1. Click the "Register" tab
2. Enter username, email, and password
3. Click "Register"
4. You'll be automatically logged in

## Step 5: Submit Your First Job

1. Click **"New Job"** in the navigation bar
2. Select read mode:
   - **Paired-end**: For paired FASTQ files (*_R1.fastq.gz and *_R2.fastq.gz)
   - **Single-end**: For single FASTQ files
3. Drag and drop your FASTQ files (or click to browse)
4. Click **"Create and Start Job"**

## Step 6: Monitor Progress

- You'll be redirected to the job monitor
- Status updates in real-time
- Watch the Nextflow log stream live
- When complete, download results from the "Results" tab

## Stopping the Service

```bash
docker-compose down
```

## Troubleshooting

### Services won't start
```bash
# Check logs
docker-compose logs backend
docker-compose logs frontend

# Restart services
docker-compose restart
```

### Can't access web interface
- Check if port 3000 is available: `lsof -i :3000`
- Try accessing: http://127.0.0.1:3000
- Check firewall settings

### Pipeline fails
- Ensure reference genome file (*.genome.fa) exists
- Check Docker has enough memory allocated
- Review job logs in the web interface

## What's Next?

- See [README_WEBSERVICE.md](README_WEBSERVICE.md) for detailed documentation
- See [TESTING.md](TESTING.md) for comprehensive testing guide
- Check API docs at http://localhost:8000/docs

## Architecture

```
┌─────────────────┐
│   Browser       │
│  (Port 3000)    │
└────────┬────────┘
         │
    ┌────▼─────┐
    │  nginx   │  Frontend
    │  React   │  (Container)
    └────┬─────┘
         │
    ┌────▼─────┐
    │ FastAPI  │  Backend
    │ SQLite   │  (Container)
    └────┬─────┘
         │
    ┌────▼─────┐
    │Nextflow  │  Pipeline
    │ Docker   │  (subprocess)
    └──────────┘
```

## Key Features

✓ Web-based interface with Material-UI
✓ User authentication and job management
✓ Drag-and-drop file upload
✓ Real-time job monitoring via WebSocket
✓ Live Nextflow log streaming
✓ Result browsing and download
✓ Supports both paired-end and single-end reads
✓ Complete RNA-seq pipeline (QC → Trim → Align → Count)

## Example Data

For testing, you can download sample RNA-seq data from:
- **NCBI SRA**: https://www.ncbi.nlm.nih.gov/sra
- **ENCODE**: https://www.encodeproject.org/
- **Use the test data generator in TESTING.md**

## Need Help?

- Full documentation: [README_WEBSERVICE.md](README_WEBSERVICE.md)
- Testing guide: [TESTING.md](TESTING.md)
- Original pipeline docs: [README.md](README.md)
- API documentation: http://localhost:8000/docs
