# RNA-seq Web Service Enhancements

This document describes the comprehensive enhancements made to transform the RNA-seq web service from a development prototype into a production-ready application.

## Phase 1: Foundation (Critical for Production)

### 1. Testing Infrastructure ✅

**Backend Tests (pytest)**
- Comprehensive test suite with 30+ tests
- Unit tests for API endpoints, authentication, and database models
- Integration tests with in-memory SQLite
- Test fixtures for users, authentication tokens, and database sessions
- Code coverage reporting with pytest-cov

**Running Tests**
```bash
cd backend
pip install -r requirements.txt
pytest
pytest --cov=app --cov-report=html
```

**Test Organization**
- `tests/test_api.py` - API endpoint tests
- `tests/test_auth.py` - Authentication and JWT tests
- `tests/test_models.py` - Database model tests
- `tests/conftest.py` - Shared fixtures and configuration

### 2. PostgreSQL Migration ✅

**Features**
- Replaced SQLite with PostgreSQL for production use
- Connection pooling (pool_size=5, max_overflow=10)
- Pre-ping for connection validation
- Connection recycling every hour

**Database Migrations (Alembic)**
- Initial schema migration (001_initial_schema.py)
- Support for both PostgreSQL and SQLite
- Automatic migration on startup

**Running Migrations**
```bash
cd backend
alembic upgrade head
alembic downgrade -1
alembic revision --autogenerate -m "description"
```

**Configuration**
```env
DATABASE_URL=postgresql://rnaseq:password@postgres:5432/rnaseq
```

### 3. Celery Job Queue System ✅

**Architecture**
- Redis as message broker and result backend
- Separate Celery worker container
- Non-blocking job submission
- Task result tracking

**Tasks Implemented**
- `run_pipeline_task` - Execute Nextflow pipeline asynchronously
- `send_email_task` - Send email notifications
- `cleanup_old_jobs_task` - Cleanup old job data

**Task Queues**
- `pipeline` - Pipeline execution tasks
- `emails` - Email notification tasks

**Configuration**
```env
CELERY_BROKER_URL=redis://redis:6379/0
CELERY_RESULT_BACKEND=redis://redis:6379/0
```

**Monitoring**
```bash
celery -A app.celery_app inspect active
celery -A app.celery_app inspect stats
```

### 4. Enhanced Security ✅

**Rate Limiting**
- Login: 5 requests/minute
- Registration: 5 requests/minute
- Job creation: 10 jobs/hour per user
- File uploads: 20 requests/minute

**CORS Configuration**
- Environment-based origin whitelist
- Explicit allowed methods and headers
- Credentials support

**Security Headers**
- X-Content-Type-Options: nosniff
- X-Frame-Options: DENY
- X-XSS-Protection: 1; mode=block
- Strict-Transport-Security: max-age=31536000
- Content-Security-Policy: Comprehensive CSP rules
- Referrer-Policy: strict-origin-when-cross-origin
- Permissions-Policy: Disabled unnecessary features

**Input Validation**
- Path validation (prevent directory traversal)
- File extension validation (.fastq.gz, .fq.gz only)
- File size limits (10GB max)
- Pydantic validators on all inputs

## Phase 2: Production Readiness

### 5. Monitoring and Logging ✅

**Structured Logging (JSON)**
- JSON-formatted logs for easy parsing
- Log levels: DEBUG, INFO, WARNING, ERROR
- Automatic timestamp and level fields
- Integration with ELK stack or similar

**Prometheus Metrics**
- HTTP request metrics (rate, duration, status codes)
- Custom application metrics
- `/metrics` endpoint for Prometheus scraping

**Enhanced Health Checks**
- Database connectivity check
- Redis connectivity check
- Disk space monitoring
- Degraded status when services are unhealthy

**Health Check Response**
```json
{
  "status": "healthy",
  "checks": {
    "database": "healthy",
    "disk_space": "healthy: 45.23GB free",
    "redis": "healthy"
  }
}
```

**Grafana Dashboard**
- Pre-configured dashboard JSON
- API request rate visualization
- Job status distribution
- Response time percentiles
- Active jobs counter

## New Services in docker-compose.yml

### PostgreSQL
```yaml
postgres:
  image: postgres:16-alpine
  environment:
    - POSTGRES_USER=rnaseq
    - POSTGRES_PASSWORD=rnaseq_dev_password
    - POSTGRES_DB=rnaseq
  volumes:
    - postgres-data:/var/lib/postgresql/data
  healthcheck:
    test: ["CMD-SHELL", "pg_isready -U rnaseq"]
```

### Redis
```yaml
redis:
  image: redis:7-alpine
  command: redis-server --appendonly yes
  volumes:
    - redis-data:/data
  healthcheck:
    test: ["CMD", "redis-cli", "ping"]
```

### Celery Worker
```yaml
celery-worker:
  build: ./backend
  command: celery -A app.celery_app worker --loglevel=info
  depends_on:
    - redis
    - postgres
```

## Environment Variables

### Required
```env
JWT_SECRET=your-strong-random-secret-key
POSTGRES_PASSWORD=your-secure-password
```

### Optional (with defaults)
```env
POSTGRES_USER=rnaseq
POSTGRES_DB=rnaseq
CORS_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
```

### Email Notifications (optional)
```env
SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USER=your-email@gmail.com
SMTP_PASSWORD=your-app-password
```

## Migration Guide

### From SQLite to PostgreSQL

1. **Backup existing data**
```bash
sqlite3 backend/rnaseq.db .dump > backup.sql
```

2. **Update environment**
```bash
cp .env.example .env
# Edit .env with your PostgreSQL credentials
```

3. **Start new services**
```bash
docker-compose up -d postgres redis
```

4. **Run migrations**
```bash
docker-compose exec backend alembic upgrade head
```

5. **Restart all services**
```bash
docker-compose restart
```

## Testing the Enhancements

### 1. Test Backend Tests
```bash
cd backend
pytest -v
```

### 2. Test PostgreSQL Connection
```bash
docker-compose exec postgres psql -U rnaseq -d rnaseq -c "SELECT version();"
```

### 3. Test Redis Connection
```bash
docker-compose exec redis redis-cli ping
```

### 4. Test Celery Workers
```bash
docker-compose exec celery-worker celery -A app.celery_app inspect ping
```

### 5. Test Health Check
```bash
curl http://localhost:8000/health | jq
```

### 6. Test Prometheus Metrics
```bash
curl http://localhost:8000/metrics
```

### 7. Test Rate Limiting
```bash
for i in {1..10}; do
  curl -X POST http://localhost:8000/api/auth/login \
    -H "Content-Type: application/json" \
    -d '{"username":"test","password":"test"}'
done
```

## Performance Improvements

### Database
- Connection pooling reduces connection overhead
- Pre-ping ensures connections are alive
- Connection recycling prevents stale connections
- Indexes on frequently queried columns

### Job Execution
- Non-blocking job submission
- Distributed task execution
- Automatic retry on failure
- Better resource utilization

### API
- Rate limiting prevents abuse
- Prometheus metrics for monitoring
- Health checks for service status
- CORS configuration reduces preflight requests

## Security Improvements

### Authentication
- JWT tokens with expiration
- Bcrypt password hashing
- Rate limiting on auth endpoints

### API Security
- Rate limiting on all endpoints
- Input validation and sanitization
- Security headers on all responses
- CORS whitelist configuration

### File Upload Security
- File extension validation
- File size limits
- Path traversal prevention
- Separate storage per job

## Monitoring Stack (Optional)

### Prometheus + Grafana Setup

Add to docker-compose.yml:
```yaml
prometheus:
  image: prom/prometheus:latest
  volumes:
    - ./monitoring/prometheus.yml:/etc/prometheus/prometheus.yml
    - prometheus-data:/prometheus
  ports:
    - "9090:9090"

grafana:
  image: grafana/grafana:latest
  volumes:
    - grafana-data:/var/lib/grafana
    - ./monitoring/grafana-dashboard.json:/etc/grafana/provisioning/dashboards/rnaseq.json
  ports:
    - "3001:3000"
  environment:
    - GF_SECURITY_ADMIN_PASSWORD=admin
```

Access:
- Prometheus: http://localhost:9090
- Grafana: http://localhost:3001 (admin/admin)

## Breaking Changes

1. **Database URL**: Now defaults to PostgreSQL instead of SQLite
2. **Job Execution**: Jobs are queued via Celery instead of direct execution
3. **CORS**: Requires explicit origin configuration (no more `allow_origins=["*"]`)
4. **Environment**: More environment variables required (see .env.example)

## Rollback Procedure

If you need to rollback to the previous version:

```bash
git revert HEAD
docker-compose down -v
docker-compose up -d
```

Note: This will lose all data in PostgreSQL. Backup first if needed.

## Future Enhancements (Phase 3+)

See the comprehensive enhancement plan for additional features:
- Email notifications
- Admin dashboard
- User quotas and limits
- Pagination and filtering
- Enhanced UI/UX
- Kubernetes deployment
- Multi-region support

## Support

For issues or questions:
1. Check logs: `docker-compose logs backend`
2. Check health: `curl http://localhost:8000/health`
3. Review this documentation
4. Check the main README_WEBSERVICE.md
