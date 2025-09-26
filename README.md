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
- For Windows users, ensure PowerShell execution policy allows script execution
- If you encounter permission issues on Windows, you may need to run PowerShell as Administrator

---

## Workflow

The following block diagram provides an overview of the **HHG-SaDAS** workflow,  
illustrating the flow of parameters, functions, and data throughout the simulation process.

<p align="center">
  <img src="workflow.svg" alt="HHG-SaDAS Workflow" width="800"/>
</p>

---


## License

This project is licensed under the [MIT License](LICENSE).
