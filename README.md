# Analysis-aware defeaturing

This repository contains the numerical illustrations published in the following papers (and some more):
- [BCV2022] "Analysis-aware defeaturing: Problem setting and _a posteriori_ estimation"\
    &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;by Annalisa Buffa, Ondine Chanon and Rafael Vázquez\
    &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;in _Mathematical Models and Methods in Applied Sciences_, 32(02), 359-402 (2022).\
    &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;Journal article: [https://doi.org/10.1142/S0218202522500099/](https://doi.org/10.1142/S0218202522500099/)\
        &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;arXiv (Open Access): [https://doi.org/10.48550/arXiv.2007.11525/](
https://doi.org/10.48550/arXiv.2007.11525/)


Libraries needed to run the examples:
  - NURBS package: [https://octave.sourceforge.io/nurbs/](https://gnu-octave.github.io/packages/nurbs/)
  - GeoPDEs 3.2.2: [http://rafavzqz.github.io/geopdes/](http://rafavzqz.github.io/geopdes/)
  


## Tests and considered geometries:

### **test01**: 
* 1 positive feature, 2D

![test01](images/test01_pos.png)


### **test02**: 
* 1 positive feature, 2D\
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^4_\varepsilon$

![test02](images/test02_pos.png)


### **test04**: 
* 1 negative feature, 2D

![test05](images/test04_neg.png)


### **test05**:
* 1 negative feature, 2D
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^2_\varepsilon$

![test05](images/test05_neg.png)


### **test06**:
* 1 negative feature and 1 positive feature, 2D
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^5_\varepsilon$

![test06](images/test06_neg_pos.png)


### **test07**:
* 1 complex feature (with a negative and a positive components), 2D
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^6_\varepsilon$

![test07](images/test07_complex.png)


### **test16**:
* 1 negative feature, 2D
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^1_\varepsilon$

![test16](images/test16_neg.png)


### **test17**:
* 1 negative feature, 2D
* [BCV2022], Section 6.1.1 where $\Omega = \Omega_s$

![test17](images/test17_neg.png)


### **test18**: 
* 1 negative feature, 2D
* [BCV2022], Section 6.1.1 where $\Omega = \Omega_c$

![test18](images/test18_neg.png)


### **test19**:
* 1 negative feature, 2D
* [BCV2022] Section 6.1.1 where $\Omega = \Omega_\star$

![test19](images/test19_neg.png)


### **test21**: 
* [BCV2022], Section 6.2.1 where $\Omega = \Omega^3_\varepsilon$
* 1 positive feature, 2D

![test21](images/test21_pos.png)


### **test30**: 
* [BCV2022], Section 6.2.2 where $\Omega = \Omega^1_\varepsilon$
* 1 negative feature, 3D

![test30](images/test30_neg.png)


### **test31**:
* [BCV2022], Section 6.2.2 where $\Omega = \Omega^3_\varepsilon$
* 1 positive feature, 3D

![test31](images/test31_pos.png)


### **test32**:
* [BCV2022], Section 6.2.2 where $\Omega = \Omega^4_\varepsilon$
* 1 positive feature, 3D

![test32](images/test32_pos.png)


### **test35**:
* [BCV2022], Section 6.2.2 where $\Omega = \Omega^2_\varepsilon$
* 1 negative feature, 3D

![test35](images/test35_neg.png)
