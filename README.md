# TFG: Post-Quantum Cryptography with Kyber and Dilithium in SageMath

## Installing SageMath on WSL (Ubuntu 22.04)

This guide explains how to install SageMath on a WSL (Windows Subsystem for Linux) environment running Ubuntu 22.04.

---

### Step 1: Update System Packages

Open your WSL terminal and run:

```bash
sudo apt update && sudo apt upgrade -y
```
### Step 2: Install SageMath

SageMath is available in the Ubuntu repositories, so installation is simple:

```bash
sudo apt install sagemath -y
```
### Step 3: Run SageMath

Once the installation is complete, you can launch SageMath by running.

We have two protocols to execute: Kyber and Dilithium

To execute the benchmarking of kyber, follow the next steps:
```bash
cd kyber
sage results.sage
```


To execute Dilithium, follow the next steps:
```bash
cd Dilithium
sage Dilithium.sage
```
