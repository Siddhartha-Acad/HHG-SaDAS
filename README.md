# HHG-SaDAS  
> For simulating Higher-order Harmonic Generation (HHG) using the Generalized Pseudospectral Method (GPSM) and the Split-Operator Method.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python](https://img.shields.io/badge/python-3.11%2B-brightgreen)
![Status](https://img.shields.io/badge/status-in%20development-orange)

---

## What is HHG-SaDAS?

**HHG-SaDAS** is a package to simulate HHG in atoms and atomic systems by solving the **3D time-dependent Schrödinger equation** using:

- **Generalized Pseudospectral Method (GPSM)**
- **Split-Operator Method**

---

This repository contains the core codes developed during **June–December 2024** as part of my *Master of Science (Research)* degree at the **Indian Institute of Technology (IIT) Mandi, India**.  
The overall tenure of my MS(R) program was **August 2022 – August 2025**.

These codes form an integral part of my MS(R) thesis and are released to ensure academic transparency and reproducibility of the research.

The original Git repository, maintained throughout the research period, is kept private.  
Whereas, this public repository is a streamlined version, providing the main useful codes in a straightforward manner.


---

### Thesis Information

- **MS(R) Thesis Title:** *Higher-Order Harmonic Generation and Harmonic-Power Enhancement in Noble-Gas Atoms Confined Inside C60*  
- **Thesis Supervisor:** Prof. Hari R. Varma  
- **Research Group:** Structure and Dynamics of Atomic Systems (SaDAS)  
- **Group Website:** [https://sadas.iitmandi.ac.in/index.php](https://sadas.iitmandi.ac.in/index.php)

---



# Installation Guide

### Prerequisites
- Git installed on your system
- Python 3.11 or higher
- PowerShell (for Windows users)


### 1. Clone the repository
```bash
git clone git@github.com:Siddhartha-Acad/HHG-SaDAS.git
cd HHG-SaDAS
```

### 2. Create and activate a virtual environment

**For Linux/macOS:**
```bash
python3 -m venv venv_HHG
source venv_HHG/bin/activate
```

**For Windows (PowerShell):**
```powershell
py -3.11 -m venv venv_HHG
.\venv_HHG\Scripts\Activate.ps1
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

## Notes
- The virtual environment step is optional but recommended to avoid package conflicts
- If you encounter permission issues on Windows, you may need to run PowerShell as Administrator

---

# Workflow

<p align="center">
  <img src="workflow.svg" alt="HHG-SaDAS Workflow" width="800"/>
</p>

### A brief description of the workflow diagram:

**Step 1**  
> Navigate to:  
> `/HHG-SaDAS/Harmonic_generation/Collocation_points/Colloc_pt_generator/`  
>  
> Here you will find the script `Algo-3_Gauss_Lobatto_colloc_pt.py`, which computes the **Gauss–Lobatto collocation points** using *Algorithm 3*, as described in Appendix A of the thesis:  
>> *“An efficient algorithm to numerically calculate the Gauss–Lobatto collocation points.”*  
>  
> After computation, the collocation points are automatically written to a `.txt` file for later use.  

**Step 2**  
> The generated collocation point data (`.txt` file) is passed into `parameters_and_functions.py`.  
> As the name suggests, this script centralizes all the required **parameters** and **Python functions** in one place.  
>  
> Users only need to modify a given parameter once inside this file. All other scripts will then read the updated values and functions directly from it.  
>  
> As indicated by the outgoing solid arrow in the workflow diagram, `parameters_and_functions.py` distributes these parameters to all subsequent scripts.  


---


## License

This project is licensed under the [MIT License](LICENSE).
