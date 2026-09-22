import axios from 'axios';

const API_BASE_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000';

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
});

api.interceptors.request.use((config) => {
  const token = localStorage.getItem('token');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  return config;
});

export const authAPI = {
  register: (username, email, password) =>
    api.post('/api/auth/register', { username, email, password }),
  
  login: (username, password) =>
    api.post('/api/auth/login', { username, password }),
};

export const jobsAPI = {
  createJob: (parameters) =>
    api.post('/api/jobs', parameters),
  
  uploadFile: (jobId, file, onProgress) => {
    const formData = new FormData();
    formData.append('file', file);
    return api.post(`/api/jobs/${jobId}/upload`, formData, {
      headers: { 'Content-Type': 'multipart/form-data' },
      onUploadProgress: onProgress,
    });
  },
  
  startJob: (jobId) =>
    api.post(`/api/jobs/${jobId}/start`),
  
  getJobs: () =>
    api.get('/api/jobs'),
  
  getJob: (jobId) =>
    api.get(`/api/jobs/${jobId}`),
  
  getJobLog: (jobId) =>
    api.get(`/api/jobs/${jobId}/log`),
  
  listResults: (jobId) =>
    api.get(`/api/jobs/${jobId}/results`),
  
  downloadResult: (jobId, filePath) =>
    `${API_BASE_URL}/api/jobs/${jobId}/results/${filePath}`,
  
  cancelJob: (jobId) =>
    api.delete(`/api/jobs/${jobId}`),
  
  getWebSocketUrl: (jobId) => {
    const wsProtocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
    const wsHost = API_BASE_URL.replace(/^https?:\/\//, '');
    return `${wsProtocol}//${wsHost}/ws/jobs/${jobId}`;
  },
};

export default api;
