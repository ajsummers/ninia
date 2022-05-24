# 0. Installing Quantum Espresso

Before getting started, you will need to install Quantum Espresso. This guide will lead you through the steps to 
install QE on an **HPC Cluster** and on **Windows**.

## 0.1 HPC Installation

For the HPC Installation, I am assuming that the system's OS is a form of Linux, and you do not have root permissions 
(this installation should work either way, but root permissions may make it simpler).

Quantum Espresso recently (as of version 7.0) started making users register before using their software. So in order to
download QE, you will have to first [register](https://www.quantum-espresso.org/register-user/). You will then be able
to [log in](https://www.quantum-espresso.org/login/) and click the first `Quantum ESPRESSO V.#.#` to download the most 
recent version of QE. You can also choose to download older versions of QE from the same page. In this walkthrough, 
V.7.0 was the most recent.

Once you download QE, move the file `qe-X.X-ReleasePack.tgz` to your desired folder in your HPC Cluster.
I decided to move it to my 'software' folder.

> ***Tip***: When you get to the download page, you can look into the page source to find the download url,
> if you can't easily transfer files between your desktop and the HPC system. In Windows, you can right click
> and select `View page source`, then press `Ctrl & F` and search for 'qe-X.X-ReleasePack.tgz'. Copy the whole url
> and use the `wget [url]` command in the HPC system to download it there. *Quantum Espresso is trying to keep track
> of their number of users, so only use this command if you have registered on their website, 
> and do not share the url.*

Next, we need to unpack the tarball:
`tar -xzvf qe-X.X-ReleasePack.tgz`

Quantum Espresso needs a few things to compile:
- A Fortran compiler compliant with F2008 standard
- A recent CMake software (v.3.14 or later)
- For parallel execution (highly recommended), an MPI-aware compiler

For these requirements, I have loaded the following modules into my HPC environment, 
which are each compatible with each other: <br>
```
module load gcc/8.4.0
module load mpich/3.2
module load cmake/3.12.4
```
> ***Notes***: My CMake version does not exactly match the requirements, but this was the latest available
> on the system I used. *It still was able to compile this version of QE.* You may need to set environmental
> variables for the modules you import. For more guidance on this, look to the Quantum Espresso user guide,
> found in the Docs on the download page.

