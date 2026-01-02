# Java Package for Heparan Sulfate Characterization

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-live-brightgreen)](https://joelrmath.github.io/heparan-sulfate/)

## Overview

This repository contains the Java implementation and R-integration workflow for the characterization of complex biochemical mixtures consisting of linear chains, e.g. Heparan Sulfate (HS). 

This work reproduces the mathematical modeling described in:
> **"Combining measurements to estimate properties and characterization extent of complex biochemical mixtures; applications to Heparan Sulfate."** > *Scientific Reports* 6, 24829 (2016). [DOI: 10.1038/srep24829](https://doi.org/10.1038/srep24829)

## Quick Links

- 🚀 **[Project Documentation & Figures](https://joelrmath.github.io/heparan-sulfate/)**: The primary landing page for model overviews and figure reproductions.
- 📚 **[API Reference (Javadoc)](https://joelrmath.github.io/heparan-sulfate/api/)**: Detailed technical documentation of the Java classes and methods.

## Getting Started

### Prerequisites
- Java JDK 11 or higher
- R (with `rJava`, `ggplot2`, and `patchwork` packages)
- Maven (for building from source)

### Installation
Clone the repository and build the JAR with dependencies:
```bash
git clone [https://github.com/JoelRMath/heparan-sulfate.git](https://github.com/JoelRMath/heparan-sulfate.git)
cd heparan-sulfate
mvn clean package
The resulting JAR will be located in the target/ directory and can be called directly from R as demonstrated in the Figure Examples.
```

Maintained by Joel R. Pradines
