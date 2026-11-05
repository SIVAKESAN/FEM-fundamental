"""FEM Learning Laboratory single-server FastAPI entry point.

Run from project root: python -m uvicorn backend.app:app --reload
"""
from pathlib import Path
from fastapi import FastAPI,HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from .models import Model
from .solver import solve_model,FemError
from .derivations import derive
from .examples import EXAMPLES

ROOT=Path(__file__).resolve().parent.parent
app=FastAPI(title='FEM Learning Laboratory',version='0.1.0',description='Educational static linear FEM for 1D/2D/3D bars, trusses, beams and frames')
app.mount('/static',StaticFiles(directory=ROOT/'static'),name='static')

@app.get('/')
def home():return FileResponse(ROOT/'static'/'index.html')

@app.get('/api/examples')
def examples():return EXAMPLES

@app.post('/api/solve')
def solve(model:Model):
    try:return solve_model(model)
    except FemError as exc:raise HTTPException(422,str(exc)) from exc

class DeriveInput(BaseModel):
    model:Model
    element_id:str

@app.post('/api/derive')
def manual(body:DeriveInput):
    try:return derive(body.model,body.element_id)
    except FemError as exc:raise HTTPException(422,str(exc)) from exc

@app.get('/api/health')
def health():return {'status':'ok'}
