#!/bin/sh
# Check if the FastAPI process is running and if the health endpoint responds
if pgrep -f "uvicorn" > /dev/null; then
  curl -s http://localhost:8000/health > /dev/null
  exit $?
else
  exit 1
fi 