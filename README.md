This repo is under active development: full installation instructions coming soon

# CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs 🐍

> 📄 **Paper**: [CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs](https://www.biorxiv.org/content/10.1101/2025.11.18.689110v1)
 
## 📝 Colab Notebook
Coming soon

## Usage

### Installation

1. **Install prerequisites**
   
   CIRPIN requires:
   
   - [Progres](https://github.com/greener-group/progres?tab=readme-ov-file)  
   - [TM-align](https://aideepmed.com/TM-align/)  
   - [Foldcomp](https://github.com/steineggerlab/foldcomp)  

3. **Updating model weights**
     
   To run CIRPIN, load ```CIRPIN/trained_models/CIRPIN_model/CIRPIN_model_5k_cp_epoch301.pt```
   
   To run Progres, load ```CIRPIN/trained_models/Progres_model/Progres_trained_model.pt```  

5. **Embedded AFDB-ClustR**
   
  The full [AFDB-ClustR](https://www.nature.com/articles/s41586-023-06510-w) embedded using CIRPIN/Progres is available for download at [link](https://huggingface.co/datasets/aidenkzj/CIRPIN/tree/main)  


## CIRPIN-DB: Accessing Databases of Circular Permutations

1. Datasets of CPs found SCOPe 40%: ```CIRPIN/scope40```
 
2. Datasets of CPs found in AFDB-ClustR: [link](https://huggingface.co/datasets/aidenkzj/CIRPIN/tree/main)


## 📄 License & Citation

**License**: MIT License - See LICENSE file for details  
**Citation**: If you use CIRPIN in your research, please cite:
```
@article {Kolodziej2025.11.18.689110,
	author = {Kolodziej, Aiden R and Abulnaga, S. Mazdak and Ovchinnikov, Sergey},
	title = {CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs},
	elocation-id = {2025.11.18.689110},
	year = {2025},
	doi = {10.1101/2025.11.18.689110},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/11/18/2025.11.18.689110},
	journal = {bioRxiv}
}
```
---

## 📧 Contact & Support

**Questions or Collaboration**: aidenkzj@mit.edu  

---
