# Getting SimNIBS GUI to work within the docker environment

1. Run the command xhost +"local:docker@" in a terminal use this for the following steps
2. Bind the container using the following specification:

	-v /usr/lib/nvidia-xxx:/usr/lib/nvidia-xxx \
	-v /usr/lib32/nvidia-xxx:/usr/lib32/nvidia-xxx \
	-e DISPLAY=unix$DISPLAY \
	--privileged \
	-v /tmp/.X11-unix:/tmp/.X11-unix

3. Within the container install pyqt5, pyopengl Python packages
4. Install qt5-default
5. export PATH="/usr/lib/nvidia-xxx/bin":$PATH
6. export LD_LIBRARY_PATH="/usr/lib/nvidia-xxx:/usr/lib32/nvidia-xxx":$LD_LIBRARY_PATH

Now you should be able to run SimNIBS GUI without any problems!

Solution from:
https://github.com/SoonminHwang/dockers/issues/1
