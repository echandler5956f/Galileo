# /home/akshay/Documents/BiQu/src/Galileo/galileo_ros/msg/EnvironmentSurface.msg

# SE3 transformation
geometry_msgs/Transform SE3_transform

# nxm dimension matrix A s.t. Ax-b <= 0 describes a point on the surface
float64[][] A

# m dimension vector b, s.t. Ax-b<=0 describes a point on the surface
float64[] b