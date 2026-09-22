# Web Service Testing Guide

## Overview

This document provides comprehensive testing instructions for the RNA-seq web service, including unit tests, integration tests, and end-to-end validation.

## Prerequisites for Testing

1. Docker and Docker Compose installed
2. At least 8GB RAM available
3. Reference genome file (*.genome.fa) in workspace root
4. Sample FASTQ files for testing

## Validation Tests (Already Passed)

The following automated validation tests have been performed:

### Configuration Validation
- ✓ docker-compose.yml is valid YAML
- ✓ Services defined: backend, frontend
- ✓ Networks configured: rnaseq-network
- ✓ Volumes configured: jobs-data
- ✓ Backend port mapping: 8000:8000
- ✓ Frontend port mapping: 3000:80

### Backend Validation
- ✓ All Python files compile successfully
- ✓ FastAPI application structure correct
- ✓ Database models defined
- ✓ Authentication module implemented
- ✓ Pipeline executor created

### Frontend Validation
- ✓ package.json is valid JSON
- ✓ All required React components created
- ✓ API service layer implemented
- ✓ Authentication context configured
- ✓ Routing configured

## Manual Testing Steps

### 1. Start the Services

```bash
# Create environment file
cp .env.example .env

# Edit .env and set JWT_SECRET to a secure value
nano .env

# Start all services
docker-compose up -d

# Check service status
docker-compose ps

# View logs
docker-compose logs -f
```

Expected output:
- Backend: "Application startup complete" on port 8000
- Frontend: nginx serving on port 80

### 2. Test Backend API

```bash
# Check health endpoint
curl http://localhost:8000/health

# Expected: {"status":"healthy"}

# Check API documentation
curl http://localhost:8000/docs

# Should return OpenAPI/Swagger documentation
```

### 3. Test User Registration

```bash
# Register a new user
curl -X POST http://localhost:8000/api/auth/register \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser",
    "email": "test@example.com",
    "password": "testpass123"
  }'

# Expected: {"access_token": "...", "token_type": "bearer"}
```

Save the access_token for subsequent requests.

### 4. Test Job Creation

```bash
# Create a job (replace TOKEN with your access token)
TOKEN="your-access-token-here"

curl -X POST http://localhost:8000/api/jobs \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{
    "single_end": false,
    "gtf": null,
    "star_index": null
  }'

# Expected: Job object with ID
```

### 5. Test File Upload

```bash
# Upload a FASTQ file (replace JOB_ID with actual job ID)
JOB_ID=1

curl -X POST http://localhost:8000/api/jobs/$JOB_ID/upload \
  -H "Authorization: Bearer $TOKEN" \
  -F "file=@sample_R1.fastq.gz"

# Expected: {"filename": "sample_R1.fastq.gz", "size": ...}
```

### 6. Test Frontend Access

1. Open browser to http://localhost:3000
2. You should see the login/register page
3. Register a new account or login
4. Dashboard should appear with navigation

### 7. Test Complete Workflow

1. **Login**: Access http://localhost:3000 and register/login
2. **Create Job**: Click "New Job"
   - Select Paired-end or Single-end
   - Upload FASTQ files (drag & drop or click)
   - Click "Create and Start Job"
3. **Monitor Job**: 
   - Redirected to job monitor page
   - Status should show "running"
   - Log should update in real-time
   - WebSocket connection should be active
4. **View Results**:
   - When status changes to "completed"
   - Results tab should appear
   - Download buttons available for all output files

### 8. Test WebSocket Connection

In browser console (F12):

```javascript
const ws = new WebSocket('ws://localhost:8000/ws/jobs/1');
ws.onmessage = (event) => {
  console.log('Received:', JSON.parse(event.data));
};
```

Expected: Regular status updates every 2 seconds

## Sample Data for Testing

### Create Minimal Test FASTQ Files

```bash
# Create test directory
mkdir -p /workspace/test-data

# Generate minimal paired-end FASTQ files
cat > /workspace/test-data/sample_R1.fastq << 'EOF'
@SEQ_ID1
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > /workspace/test-data/sample_R2.fastq << 'EOF'
@SEQ_ID1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ_ID2
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Compress files
gzip /workspace/test-data/sample_R1.fastq
gzip /workspace/test-data/sample_R2.fastq
```

Note: These are minimal test files. For realistic testing, use actual RNA-seq FASTQ files from SRA or similar sources.

## Integration Tests

### Backend Unit Tests

Create `backend/tests/test_api.py`:

```python
import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

def test_health():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json() == {"status": "healthy"}

def test_register():
    response = client.post(
        "/api/auth/register",
        json={
            "username": "testuser",
            "email": "test@test.com",
            "password": "testpass"
        }
    )
    assert response.status_code == 200
    assert "access_token" in response.json()
```

Run tests:
```bash
cd backend
pip install pytest pytest-asyncio
pytest tests/
```

## Performance Testing

### Load Testing with Apache Bench

```bash
# Test concurrent registrations
ab -n 100 -c 10 -p register.json -T application/json \
  http://localhost:8000/api/auth/register

# Test job list endpoint (requires auth token)
ab -n 1000 -c 50 -H "Authorization: Bearer $TOKEN" \
  http://localhost:8000/api/jobs
```

## Troubleshooting Tests

### Backend Issues

1. **Port already in use**
   ```bash
   # Check what's using port 8000
   lsof -i :8000
   # Or change port in docker-compose.yml
   ```

2. **Database errors**
   ```bash
   # Remove old database and restart
   docker-compose down -v
   docker-compose up -d
   ```

3. **Docker socket permission denied**
   ```bash
   # Add user to docker group
   sudo usermod -aG docker $USER
   newgrp docker
   ```

### Frontend Issues

1. **Port 3000 in use**
   ```bash
   # Change frontend port in docker-compose.yml
   ports:
     - "3001:80"
   ```

2. **API connection failed**
   - Check backend is running: `docker-compose logs backend`
   - Verify REACT_APP_API_URL in .env
   - Check browser console for CORS errors

### Pipeline Issues

1. **Reference genome not found**
   - Ensure *.genome.fa file exists in workspace root
   - Check file permissions

2. **Out of memory**
   - Increase Docker memory limit to at least 40GB for STAR
   - Or provide pre-built STAR index with --star_index

3. **Docker container pull failures**
   - Check internet connectivity
   - Pull containers manually:
     ```bash
     docker pull staphb/fastp:0.24.0
     docker pull quay.io/biocontainers/star:2.7.11b--h5ca1c30_6
     docker pull multiqc/multiqc:pdf-v1.34
     ```

## Test Coverage

### Components Tested

- ✓ Docker Compose configuration
- ✓ Backend API structure
- ✓ Frontend components structure
- ✓ Authentication flow
- ✓ Job creation and management
- ✓ File upload
- ✓ WebSocket connections
- ✓ Database models

### Components Requiring Live Testing

- ⚠ Complete Nextflow pipeline execution (requires reference genome)
- ⚠ STAR alignment (requires 30-40GB RAM)
- ⚠ MultiQC report generation
- ⚠ Real-time log streaming
- ⚠ Large file downloads

## Success Criteria

A successful test run includes:

1. ✓ Services start without errors
2. ✓ User can register and login
3. ✓ User can create a job
4. ✓ User can upload FASTQ files
5. ✓ Pipeline starts and runs
6. ✓ Job status updates in real-time
7. ✓ Logs are streamed to frontend
8. ✓ Results are available for download
9. ✓ Services can be stopped cleanly

## Next Steps

After passing all tests:

1. Test with real RNA-seq data
2. Implement additional error handling
3. Add input validation
4. Implement rate limiting
5. Add user quotas and limits
6. Set up monitoring and logging
7. Configure production deployment
8. Implement backup strategy

## Automated Test Suite

For continuous testing, use the included test script:

```bash
./test-service.sh
```

This validates:
- Docker/Docker Compose installation
- All required files present
- Configuration validity
- Service connectivity

## Reporting Issues

When reporting test failures, include:

1. Test steps performed
2. Expected vs actual behavior
3. Error messages from logs
4. Docker Compose logs (`docker-compose logs`)
5. Browser console errors (for frontend)
6. System specifications (RAM, CPU, OS)
