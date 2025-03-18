

# Overall storyline
## NLP
  
* [x] Scaling with regards to elements
    -  With full, half and no initialization
* [x] Scaling with regards to datapoints
    -  With full, half and no initialization
* [x] Computational cost of extending to nonlinear strain
* [x] L1 vs L2 strain

## Alternating dirs

* [ ] Global optima comparison with NLP - with noise
    - [ ] Linear strain
        - [ ] With and without initialization
    - [ ] Nonlinear strain
        - [ ] With and without initialization
* [ ] Scaling with regards to elements
    - [ ] Without initialization
    - [ ] Linear and nonlinear strain
* [ ] Scaling with regards to datapoints
    - [ ] With and without initialization
    - [ ] Linear and nonlinear strain
    

# Questions

* Is "use_data_bounds" of any interest, or is it better to just explain it?

# TODO

## LP solver

* Update tests to include both types of initialization
* Update tests to check that they all get the same results
* Write a good test to check whether different E, S points are choosen
    - Got them to choose different points when increasing to 300 points with noise, without noise, no difference up to 400 points
* Measure computational cost of extending to nonlinear strain

## Alternating direction solver

* Compare with and without initialization
    - Check cost
    - Compare iterations

* Compare 