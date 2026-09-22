import React, { useState, useEffect, useRef } from 'react';
import {
  Paper,
  Typography,
  Box,
  Chip,
  Button,
  Card,
  CardContent,
  Grid,
  Alert,
  CircularProgress,
  Tabs,
  Tab,
} from '@mui/material';
import { Refresh, Download, Cancel } from '@mui/icons-material';
import { useParams, useNavigate } from 'react-router-dom';
import { jobsAPI } from '../services/api';

function JobMonitor() {
  const { jobId } = useParams();
  const [job, setJob] = useState(null);
  const [log, setLog] = useState('');
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(true);
  const [tab, setTab] = useState(0);
  const wsRef = useRef(null);
  const navigate = useNavigate();

  useEffect(() => {
    loadJob();
    loadLog();

    const ws = new WebSocket(jobsAPI.getWebSocketUrl(jobId));
    
    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      setJob((prev) => ({ ...prev, status: data.status }));
      if (data.log) {
        setLog(data.log);
      }
    };

    ws.onerror = (error) => {
      console.error('WebSocket error:', error);
    };

    wsRef.current = ws;

    return () => {
      if (wsRef.current) {
        wsRef.current.close();
      }
    };
  }, [jobId]);

  useEffect(() => {
    if (job && job.status === 'completed') {
      loadResults();
    }
  }, [job]);

  const loadJob = async () => {
    try {
      const response = await jobsAPI.getJob(jobId);
      setJob(response.data);
    } catch (error) {
      console.error('Failed to load job:', error);
    } finally {
      setLoading(false);
    }
  };

  const loadLog = async () => {
    try {
      const response = await jobsAPI.getJobLog(jobId);
      setLog(response.data.log);
    } catch (error) {
      console.error('Failed to load log:', error);
    }
  };

  const loadResults = async () => {
    try {
      const response = await jobsAPI.listResults(jobId);
      setResults(response.data.files);
    } catch (error) {
      console.error('Failed to load results:', error);
    }
  };

  const handleCancel = async () => {
    try {
      await jobsAPI.cancelJob(jobId);
      loadJob();
    } catch (error) {
      console.error('Failed to cancel job:', error);
    }
  };

  const getStatusColor = (status) => {
    const colors = {
      pending: 'default',
      running: 'primary',
      completed: 'success',
      failed: 'error',
      cancelled: 'warning',
    };
    return colors[status] || 'default';
  };

  if (loading) {
    return (
      <Box sx={{ display: 'flex', justifyContent: 'center', mt: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  if (!job) {
    return <Alert severity="error">Job not found</Alert>;
  }

  return (
    <Box>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 3 }}>
        <Typography variant="h5">Job #{job.id}</Typography>
        <Box>
          <Button
            startIcon={<Refresh />}
            onClick={loadJob}
            sx={{ mr: 1 }}
          >
            Refresh
          </Button>
          {job.status === 'running' && (
            <Button
              startIcon={<Cancel />}
              onClick={handleCancel}
              color="error"
            >
              Cancel
            </Button>
          )}
          <Button onClick={() => navigate('/dashboard')} sx={{ ml: 1 }}>
            Back to Jobs
          </Button>
        </Box>
      </Box>

      <Grid container spacing={3} sx={{ mb: 3 }}>
        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" gutterBottom>
                Status
              </Typography>
              <Chip
                label={job.status}
                color={getStatusColor(job.status)}
              />
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" gutterBottom>
                Mode
              </Typography>
              <Typography variant="h6">
                {job.parameters.single_end ? 'Single-end' : 'Paired-end'}
              </Typography>
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" gutterBottom>
                Created
              </Typography>
              <Typography variant="body2">
                {new Date(job.created_at).toLocaleString()}
              </Typography>
            </CardContent>
          </Card>
        </Grid>
        <Grid item xs={12} md={3}>
          <Card>
            <CardContent>
              <Typography color="text.secondary" gutterBottom>
                Updated
              </Typography>
              <Typography variant="body2">
                {new Date(job.updated_at).toLocaleString()}
              </Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      {job.error_message && (
        <Alert severity="error" sx={{ mb: 3 }}>
          {job.error_message}
        </Alert>
      )}

      <Paper sx={{ mb: 3 }}>
        <Tabs value={tab} onChange={(e, v) => setTab(v)}>
          <Tab label="Log" />
          {job.status === 'completed' && <Tab label="Results" />}
        </Tabs>

        <Box sx={{ p: 2 }}>
          {tab === 0 && (
            <Box
              sx={{
                bgcolor: 'grey.900',
                color: 'grey.100',
                p: 2,
                borderRadius: 1,
                fontFamily: 'monospace',
                fontSize: '0.875rem',
                maxHeight: '500px',
                overflow: 'auto',
                whiteSpace: 'pre-wrap',
              }}
            >
              {log || 'No log available yet...'}
            </Box>
          )}

          {tab === 1 && job.status === 'completed' && (
            <Box>
              {results.length === 0 ? (
                <Typography>Loading results...</Typography>
              ) : (
                <Box>
                  <Typography variant="h6" gutterBottom>
                    Results Files ({results.length})
                  </Typography>
                  {results.map((file, index) => (
                    <Box
                      key={index}
                      sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        p: 1,
                        borderBottom: '1px solid',
                        borderColor: 'divider',
                      }}
                    >
                      <Box>
                        <Typography variant="body2">{file.path}</Typography>
                        <Typography variant="caption" color="text.secondary">
                          {(file.size / 1024 / 1024).toFixed(2)} MB
                        </Typography>
                      </Box>
                      <Button
                        startIcon={<Download />}
                        size="small"
                        href={jobsAPI.downloadResult(jobId, file.path)}
                        target="_blank"
                      >
                        Download
                      </Button>
                    </Box>
                  ))}
                </Box>
              )}
            </Box>
          )}
        </Box>
      </Paper>
    </Box>
  );
}

export default JobMonitor;
