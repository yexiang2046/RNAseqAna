import React, { useState } from 'react';
import {
  Paper,
  Typography,
  Box,
  TextField,
  FormControl,
  FormLabel,
  RadioGroup,
  FormControlLabel,
  Radio,
  Button,
  List,
  ListItem,
  ListItemText,
  IconButton,
  Alert,
  CircularProgress,
  LinearProgress,
} from '@mui/material';
import { Delete, CloudUpload } from '@mui/icons-material';
import { useDropzone } from 'react-dropzone';
import { useNavigate } from 'react-router-dom';
import { jobsAPI } from '../services/api';

function JobUpload() {
  const [mode, setMode] = useState('paired');
  const [gtf, setGtf] = useState('');
  const [starIndex, setStarIndex] = useState('');
  const [files, setFiles] = useState([]);
  const [uploading, setUploading] = useState(false);
  const [uploadProgress, setUploadProgress] = useState({});
  const [error, setError] = useState('');
  const navigate = useNavigate();

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    accept: {
      'application/gzip': ['.gz'],
      'application/x-gzip': ['.fastq.gz'],
    },
    onDrop: (acceptedFiles) => {
      setFiles([...files, ...acceptedFiles]);
    },
  });

  const removeFile = (index) => {
    setFiles(files.filter((_, i) => i !== index));
  };

  const handleSubmit = async () => {
    if (files.length === 0) {
      setError('Please upload at least one FASTQ file');
      return;
    }

    setUploading(true);
    setError('');

    try {
      const jobResponse = await jobsAPI.createJob({
        single_end: mode === 'single',
        gtf: gtf || null,
        star_index: starIndex || null,
      });

      const jobId = jobResponse.data.id;

      for (let i = 0; i < files.length; i++) {
        const file = files[i];
        await jobsAPI.uploadFile(jobId, file, (progressEvent) => {
          const percentCompleted = Math.round(
            (progressEvent.loaded * 100) / progressEvent.total
          );
          setUploadProgress((prev) => ({
            ...prev,
            [file.name]: percentCompleted,
          }));
        });
      }

      await jobsAPI.startJob(jobId);

      navigate(`/dashboard/jobs/${jobId}`);
    } catch (error) {
      setError(error.response?.data?.detail || 'Failed to create job');
      setUploading(false);
    }
  };

  return (
    <Paper sx={{ p: 3 }}>
      <Typography variant="h5" gutterBottom>
        Create New Job
      </Typography>

      {error && <Alert severity="error" sx={{ mb: 2 }}>{error}</Alert>}

      <Box sx={{ mb: 3 }}>
        <FormControl component="fieldset">
          <FormLabel component="legend">Read Mode</FormLabel>
          <RadioGroup
            row
            value={mode}
            onChange={(e) => setMode(e.target.value)}
          >
            <FormControlLabel value="paired" control={<Radio />} label="Paired-end" />
            <FormControlLabel value="single" control={<Radio />} label="Single-end" />
          </RadioGroup>
        </FormControl>
      </Box>

      <TextField
        fullWidth
        label="GTF File Path (optional)"
        value={gtf}
        onChange={(e) => setGtf(e.target.value)}
        sx={{ mb: 2 }}
        helperText="Leave empty to use default"
      />

      <TextField
        fullWidth
        label="STAR Index Path (optional)"
        value={starIndex}
        onChange={(e) => setStarIndex(e.target.value)}
        sx={{ mb: 3 }}
        helperText="Leave empty to build index automatically"
      />

      <Box
        {...getRootProps()}
        sx={{
          border: '2px dashed',
          borderColor: isDragActive ? 'primary.main' : 'grey.400',
          borderRadius: 2,
          p: 3,
          textAlign: 'center',
          cursor: 'pointer',
          bgcolor: isDragActive ? 'action.hover' : 'background.paper',
          mb: 2,
        }}
      >
        <input {...getInputProps()} />
        <CloudUpload sx={{ fontSize: 48, color: 'action.active', mb: 1 }} />
        <Typography>
          {isDragActive
            ? 'Drop FASTQ files here'
            : 'Drag and drop FASTQ files here, or click to select'}
        </Typography>
        <Typography variant="caption" color="text.secondary">
          Supports .fastq.gz files
        </Typography>
      </Box>

      {files.length > 0 && (
        <Box sx={{ mb: 2 }}>
          <Typography variant="subtitle2" gutterBottom>
            Files to upload ({files.length})
          </Typography>
          <List>
            {files.map((file, index) => (
              <ListItem
                key={index}
                secondaryAction={
                  !uploading && (
                    <IconButton edge="end" onClick={() => removeFile(index)}>
                      <Delete />
                    </IconButton>
                  )
                }
              >
                <ListItemText
                  primary={file.name}
                  secondary={`${(file.size / 1024 / 1024).toFixed(2)} MB`}
                />
                {uploading && uploadProgress[file.name] !== undefined && (
                  <Box sx={{ width: '100px', ml: 2 }}>
                    <LinearProgress
                      variant="determinate"
                      value={uploadProgress[file.name]}
                    />
                  </Box>
                )}
              </ListItem>
            ))}
          </List>
        </Box>
      )}

      <Button
        fullWidth
        variant="contained"
        size="large"
        onClick={handleSubmit}
        disabled={uploading || files.length === 0}
        startIcon={uploading && <CircularProgress size={20} />}
      >
        {uploading ? 'Creating Job...' : 'Create and Start Job'}
      </Button>
    </Paper>
  );
}

export default JobUpload;
