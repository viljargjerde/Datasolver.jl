

# Overall storyline
## NLP
  
* [x] Scaling with regards to elements
    -  With full, half and no initialization
* [x] Scaling with regards to datapoints
    -  With full, half and no initialization
* [x] Computational cost of extending to nonlinear strain
* [x] L1 vs L2 strain

## Alternating dirs

* [x] Global optima comparison with NLP - without noise
    - Linear strain
        - With and without initialization
    - Nonlinear strain
        - With and without initialization
* [x] Global optima comparison with NLP - with noise
    - Linear strain
        - With and without initialization
    - Nonlinear strain
        - With and without initialization
* [ ] Scaling with regards to elements
    - [ ] Without initialization
    - [ ] Linear and nonlinear strain
* [ ] Scaling with regards to datapoints
    - [ ] With and without initialization
    - [ ] Linear and nonlinear strain
    

# Questions

* Is "use_data_bounds" of any interest, or is it better to just explain it?


# Notes

* It is necessary to use different units to achieve numerical stability for the LP solver.