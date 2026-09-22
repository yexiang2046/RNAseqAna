import pytest
from fastapi import status


class TestHealthEndpoint:
    def test_health_check(self, client):
        response = client.get("/health")
        assert response.status_code == status.HTTP_200_OK
        assert response.json() == {"status": "healthy"}

    def test_root_endpoint(self, client):
        response = client.get("/")
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert "message" in data
        assert "version" in data


class TestAuthenticationEndpoints:
    def test_register_new_user(self, client):
        response = client.post(
            "/api/auth/register",
            json={
                "username": "newuser",
                "email": "newuser@example.com",
                "password": "password123"
            }
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert "access_token" in data
        assert data["token_type"] == "bearer"

    def test_register_duplicate_username(self, client, test_user):
        response = client.post(
            "/api/auth/register",
            json={
                "username": "testuser",
                "email": "different@example.com",
                "password": "password123"
            }
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST
        assert "already registered" in response.json()["detail"].lower()

    def test_register_duplicate_email(self, client, test_user):
        response = client.post(
            "/api/auth/register",
            json={
                "username": "differentuser",
                "email": "test@example.com",
                "password": "password123"
            }
        )
        assert response.status_code == status.HTTP_400_BAD_REQUEST

    def test_login_success(self, client, test_user):
        response = client.post(
            "/api/auth/login",
            json={"username": "testuser", "password": "testpass123"}
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert "access_token" in data
        assert data["token_type"] == "bearer"

    def test_login_wrong_password(self, client, test_user):
        response = client.post(
            "/api/auth/login",
            json={"username": "testuser", "password": "wrongpassword"}
        )
        assert response.status_code == status.HTTP_401_UNAUTHORIZED

    def test_login_nonexistent_user(self, client):
        response = client.post(
            "/api/auth/login",
            json={"username": "nonexistent", "password": "password"}
        )
        assert response.status_code == status.HTTP_401_UNAUTHORIZED


class TestJobEndpoints:
    def test_create_job(self, client, auth_headers):
        response = client.post(
            "/api/jobs",
            json={"single_end": False, "gtf": None, "star_index": None},
            headers=auth_headers
        )
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["status"] == "pending"
        assert data["parameters"]["single_end"] is False

    def test_create_job_unauthorized(self, client):
        response = client.post(
            "/api/jobs",
            json={"single_end": False}
        )
        assert response.status_code == status.HTTP_403_FORBIDDEN

    def test_list_jobs_empty(self, client, auth_headers):
        response = client.get("/api/jobs", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK
        assert response.json() == []

    def test_list_jobs_with_data(self, client, auth_headers):
        client.post(
            "/api/jobs",
            json={"single_end": True},
            headers=auth_headers
        )
        
        response = client.get("/api/jobs", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK
        jobs = response.json()
        assert len(jobs) == 1
        assert jobs[0]["status"] == "pending"

    def test_get_job_by_id(self, client, auth_headers):
        create_response = client.post(
            "/api/jobs",
            json={"single_end": False},
            headers=auth_headers
        )
        job_id = create_response.json()["id"]

        response = client.get(f"/api/jobs/{job_id}", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK
        data = response.json()
        assert data["id"] == job_id
        assert data["status"] == "pending"

    def test_get_nonexistent_job(self, client, auth_headers):
        response = client.get("/api/jobs/99999", headers=auth_headers)
        assert response.status_code == status.HTTP_404_NOT_FOUND

    def test_get_other_users_job(self, client, auth_headers, db_session):
        from app import models
        from app.auth import get_password_hash
        
        other_user = models.User(
            username="otheruser",
            email="other@example.com",
            hashed_password=get_password_hash("pass")
        )
        db_session.add(other_user)
        db_session.commit()
        
        other_job = models.Job(
            user_id=other_user.id,
            status="pending",
            parameters={"single_end": False}
        )
        db_session.add(other_job)
        db_session.commit()

        response = client.get(f"/api/jobs/{other_job.id}", headers=auth_headers)
        assert response.status_code == status.HTTP_404_NOT_FOUND

    def test_cancel_job(self, client, auth_headers):
        create_response = client.post(
            "/api/jobs",
            json={"single_end": False},
            headers=auth_headers
        )
        job_id = create_response.json()["id"]

        response = client.delete(f"/api/jobs/{job_id}", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK

    def test_get_job_log(self, client, auth_headers):
        create_response = client.post(
            "/api/jobs",
            json={"single_end": False},
            headers=auth_headers
        )
        job_id = create_response.json()["id"]

        response = client.get(f"/api/jobs/{job_id}/log", headers=auth_headers)
        assert response.status_code == status.HTTP_200_OK
        assert "log" in response.json()
