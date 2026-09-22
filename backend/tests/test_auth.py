import pytest
from app.auth import (
    verify_password,
    get_password_hash,
    create_access_token,
    decode_token
)
from fastapi import HTTPException


class TestPasswordHashing:
    def test_hash_password(self):
        password = "mySecurePassword123"
        hashed = get_password_hash(password)
        assert hashed != password
        assert len(hashed) > 0

    def test_verify_correct_password(self):
        password = "testPassword"
        hashed = get_password_hash(password)
        assert verify_password(password, hashed) is True

    def test_verify_wrong_password(self):
        password = "testPassword"
        hashed = get_password_hash(password)
        assert verify_password("wrongPassword", hashed) is False


class TestJWTTokens:
    def test_create_token(self):
        data = {"sub": "testuser"}
        token = create_access_token(data)
        assert isinstance(token, str)
        assert len(token) > 0

    def test_decode_valid_token(self):
        username = "testuser"
        token = create_access_token({"sub": username})
        decoded_username = decode_token(token)
        assert decoded_username == username

    def test_decode_invalid_token(self):
        with pytest.raises(HTTPException) as exc_info:
            decode_token("invalid.token.here")
        assert exc_info.value.status_code == 401

    def test_token_without_sub(self):
        token = create_access_token({"other": "data"})
        with pytest.raises(HTTPException) as exc_info:
            decode_token(token)
        assert exc_info.value.status_code == 401
