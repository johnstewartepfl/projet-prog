# SUPER PROJET PROG
### chemical equation equilibrer 

programation projet with Fane Shala, John Stewart and Albéric Vigne 

- [ ] Takes reactants and products as inpts 
- [ ] Atom check
- [ ] Stoechiometry calculation
- [ ] Display the reaction
- [ ] Calculate the standard free gibbs energy 

usefull links for better understanding: 
 * [hello](coucou.fr)


![Project Logo](projet-prog.png)

![Coverage Status](https://raw.githubusercontent.com/pschwllr/minimal_project/main/assets/coverage-badge.svg)


<br>


Chemical equation balancer

## Goal

> Balance cehmical equation, gives a representation and the standard free gibbs energy (TODO show in a very small amount of space the **MOST** useful thing your package can do.)


## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n ch200 python=3.10 
```

```
conda activate ch200
```

If you need jupyter lab, install it 

```
(ch200) $ pip install jupyterlab
```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:johnstewartepfl/projet-prog'.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:pschwllr/ch200.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(ch200) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```

### Generate coverage badge

Works after running `tox`

```
(conda_env) $ pip install "genbadge[coverage]"
(conda_env) $ genbadge coverage -i coverage.xml
```



