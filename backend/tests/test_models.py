import pytest
from datetime import datetime
from app import models
from app.auth import get_password_hash


class TestUserModel:
    def test_create_user(self, db_session):
        user = models.User(
            username="testuser",
            email="test@example.com",
            hashed_password=get_password_hash("password")
        )
        db_session.add(user)
        db_session.commit()
        db_session.refresh(user)

        assert user.id is not None
        assert user.username == "testuser"
        assert user.email == "test@example.com"
        assert user.created_at is not None
        assert isinstance(user.created_at, datetime)

    def test_user_unique_username(self, db_session):
        user1 = models.User(
            username="same",
            email="email1@example.com",
            hashed_password="hash1"
        )
        user2 = models.User(
            username="same",
            email="email2@example.com",
            hashed_password="hash2"
        )
        
        db_session.add(user1)
        db_session.commit()
        
        db_session.add(user2)
        with pytest.raises(Exception):
            db_session.commit()

    def test_user_unique_email(self, db_session):
        user1 = models.User(
            username="user1",
            email="same@example.com",
            hashed_password="hash1"
        )
        user2 = models.User(
            username="user2",
            email="same@example.com",
            hashed_password="hash2"
        )
        
        db_session.add(user1)
        db_session.commit()
        
        db_session.add(user2)
        with pytest.raises(Exception):
            db_session.commit()


class TestJobModel:
    def test_create_job(self, db_session, test_user):
        job = models.Job(
            user_id=test_user.id,
            status="pending",
            parameters={"single_end": False, "gtf": None}
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)

        assert job.id is not None
        assert job.user_id == test_user.id
        assert job.status == "pending"
        assert job.parameters["single_end"] is False
        assert job.created_at is not None
        assert job.updated_at is not None

    def test_job_relationship_with_user(self, db_session, test_user):
        job = models.Job(
            user_id=test_user.id,
            status="running",
            parameters={}
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)

        assert job.user.username == test_user.username
        assert job in test_user.jobs

    def test_job_status_values(self, db_session, test_user):
        statuses = ["pending", "running", "completed", "failed", "cancelled"]
        
        for status_val in statuses:
            job = models.Job(
                user_id=test_user.id,
                status=status_val,
                parameters={}
            )
            db_session.add(job)
            db_session.commit()
            db_session.refresh(job)
            assert job.status == status_val

    def test_job_with_results_path(self, db_session, test_user):
        job = models.Job(
            user_id=test_user.id,
            status="completed",
            parameters={},
            results_path="/workspace/jobs/1/results"
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)

        assert job.results_path == "/workspace/jobs/1/results"

    def test_job_with_error_message(self, db_session, test_user):
        error_msg = "Pipeline failed: out of memory"
        job = models.Job(
            user_id=test_user.id,
            status="failed",
            parameters={},
            error_message=error_msg
        )
        db_session.add(job)
        db_session.commit()
        db_session.refresh(job)

        assert job.error_message == error_msg

    def test_job_cascade_delete_user(self, db_session):
        user = models.User(
            username="deleteuser",
            email="delete@example.com",
            hashed_password="hash"
        )
        db_session.add(user)
        db_session.commit()

        job = models.Job(
            user_id=user.id,
            status="pending",
            parameters={}
        )
        db_session.add(job)
        db_session.commit()

        job_id = job.id
        
        db_session.delete(user)
        db_session.commit()

        deleted_job = db_session.query(models.Job).filter(
            models.Job.id == job_id
        ).first()
        
        assert deleted_job is None
