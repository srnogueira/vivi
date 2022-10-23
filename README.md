# Vivi

Vivi is a framework for energy integration optimization. It contains structures and functions to facilitate the selection/sizing of technologies (process and utilities), under constraints of mass, energy and heat cascade (pinch analysis) balances.

## Optimization problem
Vivi will write a linear optimization problem as:

$$ \max_{\gamma^\omega,in_i,out_i,R_k} \left [ \sum_o^{O} out_o \cdot v_o - \sum_i^{I} in_i \cdot v_i \right ] $$

$s.t.$

$$ in_i + \sum_{\omega}^{\Omega} \gamma^\omega\left (T_{out,i}^{\omega} - T_{in,i}^{\omega} \right) = out_i \text{ } \forall \text{ } i \in \text{Resources} $$

$$ R_{k-1} + \sum_{\omega}^{\Omega} \left (\gamma^{\omega} \sum_n^N \dot Q_{n,k}^{\omega} \right) = R_{k} \text{ } \forall \text{ } k \in \text{ Heat cascade}$$

In which,

$$ \gamma^\omega \geq 0 $$

$$ in_i  \begin{cases} 
            \geq 0 &\text{if } in_i = Inf \\
            = in_i &\text{otherwise}
        \end{cases} $$


$$ out_i  \begin{cases} 
            \geq 0 &\text{if } out_i = Inf \\
            = out_i &\text{otherwise}
        \end{cases} $$

$$ R_k  \begin{cases} 
            = 0 &\text{if } k = 0 \text{ or } k = K \\
            \geq 0 &\text{otherwise}
        \end{cases} $$

| Name | Description |
| -------- | ----------- |
| $in_i$ and $out_i$ | amount of resource input and output "i", respectively | 
| $v_i$ | specific value of resource "i" |
| $\gamma^\omega$ | size factor for technology $\omega$ |
| $T_{in,i}^{\omega}$ and $T_{out,i}^{\omega}$ | input and output amount of resource "i" for technology $\omega$, respectively |
| $R_k$ | Net heat of temperature interval "k" in Heat cascade |
| $\dot Q_{n,k}^\omega $ | Heat transfer of stream "n" in temperature interval "k" for technology $\omega$ |
| $i$ and $I$ | input number and total number of inputs, respectively |
| $o$ and $O$ | output number and total number of outputs, respectively |
| $n$ and $N$ | stream number and total number of streams, respectively |
| $k$ and $K$ | temperature interval number and total number of temperature intervals, respectively |

In simple terms, vivi will determine the set of size factors $(\gamma^\omega)$, inputs $(in_i)$ and outputs $(out_o)$ that optimize the value of outputs, discounting the inputs value, and respecting the balance of resources and heat cascade (first and second constraint, respectively). If the amount of a certain input or output is fixed (e.g., $in_i\neq Inf$), this information it is also accounted as an additional constraint.

**Note**: "value" can be a monetary figure, energy, exergy, carbon, amoung others.

## How to use?
See the jupyter notebook "tutorial.ipynb" in the folder "examples"

## Citation
Nakashima, R. N. (2022). Modelling, simulation and optimization of biogas conversion routes integrated with fuel cell technology. Doctoral Thesis, Escola Politécnica, University of São Paulo, São Paulo. https://doi.org/10.11606/T.3.2022.tde-26082022-081436


### Bibtex

    
    @phdthesis{nakashima_modelling_2022,
	    address = {São Paulo},
	    title = {Modelling, simulation and optimization of biogas conversion routes integrated with fuel cell technology.},
	    url = {https://www.teses.usp.br/teses/disponiveis/3/3150/tde-26082022-081436/},
	    language = {en},
	    urldate = {2022-09-05},
	    school = {Universidade de São Paulo},
	    author = {Nakashima, Rafael Nogueira},
	    month = feb,
	    year = {2022},
	    doi = {10.11606/T.3.2022.tde-26082022-081436},
    }


## Contact
For doubts, comments and requests:

> Rafael Nogueira Nakashima (rafnn@dtu.dk)

## References

Papoulias, S. A., & Grossmann, I. E. (1983). A structural optimization approach in process synthesis—II. Computers & Chemical Engineering, 7(6), 707–721. https://doi.org/10.1016/0098-1354(83)85023-6

Marechal, F., & Kalitventzeff, B. (1996). Targeting the minimum cost of energy requirements: A new graphical technique for evaluating the integration of utility systems. Computers & Chemical Engineering, 20, S225–S230. https://doi.org/10.1016/0098-1354(96)00048-8

Maréchal, F., & Kalitventzeff, B. (1997). Effect modelling and optimization, a new methodology for combined energy and environment synthesis of industrial processes. Applied Thermal Engineering, 17(8–10), 981–992. https://doi.org/10.1016/S1359-4311(96)00079-8

Kemp, I. C. (2006). Pinch Analysis and Process Integration. Elsevier. https://doi.org/10.1016/B978-0-7506-8260-2.X5001-9

