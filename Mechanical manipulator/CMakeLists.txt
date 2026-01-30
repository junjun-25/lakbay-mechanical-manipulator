# For more information about build system see
# https://docs.espressif.com/projects/esp-idf/en/latest/api-guides/build-system.html
# The following five lines of boilerplate have to be in your project's
# CMakeLists in this exact order for cmake to work correctly
cmake_minimum_required(VERSION 3.5)

set(EXTRA_COMPONENT_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/components/micro_ros_espidf_component
)


include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(Robot_Arm)

## sudo docker run -it --rm -v /dev:/dev --privileged --net=host microros/micro-ros-agent:jazzy serial --dev /dev/ttyUSB0 -b 921600 -v6
## sudo docker run -it --rm -v /dev:/dev --privileged --net=host microros/micro-ros-agent:jazzy serial --dev /dev/ttyACM0 -v6

#ros2 topic pub -r 1 /traj_point trajectory_msgs/msg/JointTrajectoryPoint "{positions: [0.0, 0.5, 1.0, -0.5, 0.25, -1.0], velocities: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], accelerations: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], effort: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], time_from_start: {sec: 0, nanosec: 10000000}}"

#ros2 topic pub -r 1 /robot_pose geometry_msgs/msg/Pose "position: x: 0.5 y: 0.0 z: 0.2 orientation: x: 0.0 y: 0.0 z: 0.0 w: 1.0"

#ros2 topic pub -r 1 /traj_point std_msgs/msg/Float32MultiArray "data: [0.0, 0.5, 1.0, -0.5, 0.25, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]"


## ros2 topic pub -r 50 /target_pose geometry_msgs/msg/Pose "{position: {x: 1300.0, y: 0.0, z: 100.0}, orientation: {x: 0.0, y: 0.0, z: 0.0, w: 1.0}}"
