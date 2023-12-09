# Original Galileo Python Tests

These were the original python experiments which were used as a proof of concept for Galileo. These tests were inspired by the JNRH'2023 Pinocchio Tutorial which showcased how to combine Casadi and Pinocchio, but the project quickly grew far beyond the scope of what was presented in these tests. We now have a higher focus on MPC and contact decisions in the main Galileo framework.

## Getting the Python Tests to Work

Start by cloning this repository:
```
git clone https://github.com/echandler5956f/Galileo.git
cd Galileo
```

After this, you need to install the relevant dependencies. Note that the python tests require Pinocchio V3, since it has the python bindings with Casadi.

### pip

1. create an environment:

    `python3 -m venv .venv`

2. activate it:

    `source .venv/bin/activate`

3. update pip:

    `pip install -U pip`

4. install dependencies:

    `pip install casadi pin3x-jnrh2023 frozendict meshcat scipy matplotlib`

Now, you can run the various tests

```bash
python3 python/huron_centroidal_v4.py
```
or
```bash
python3 python/test_huron_centroidal.py
```

Note: Since I have access to certain HSL solvers, I am using `ma97` for many of the tests. For this public release, I have this solver selection commented out, but if you have access to `ma97`, I highly recommend uncommenting it, because it substantially improves performance.

## TODO: 

Add explanation for how to use the whole body controller with HURON in `python/test_quad_prog.py`