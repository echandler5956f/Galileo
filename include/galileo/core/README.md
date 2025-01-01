# galileo/core
The `galileo/core` folder contains the CRTP base classes necessary to define and solve optimal control problems.

The base classes are broken up into `...ModelBase` and `...DataBase` classes. 

`...Model...` classes hold information that is quasi-static while solving an optimal control problem. For instance, a `SegmentModelERK` should not change its underlying integration model (which is defined by a Butcher tableu) over the course of an optimal control problem—but between consecutive optimal control problem instances, it may change (for instance, we may initially want to find a course-grained solution with RK2, and then want a more accurate solution with RK4). 

On the other hand, `...Data...` classes hold information that is expected to change or be updated over the course of solving an optimal control problem. For instance, a `NodeData...` class stores derivatives of the constraints at a certain node calculated from `NodeModel...::calcDiff(x, u)`.

`galileo/core` also contains `visitors` used for iterating over heterogeneous sets of derived classes from a certain `...Base` CRTP class. The heterogeneous sets are defined using `boost::variant`. Unfortunately, due to the way variants work, we need to declare all possible types that the variant might encounter ahead of time.