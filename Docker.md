Examples of special usages
----------------------
- [Install Docker](#install-docker)
- [Create Docker image for grabb](#create-docker-image-for-grabb)
- [Run the Docker image and attach local folders](#run-the-docker-image-and-attach-local-folders)
- [Run the Docker image detached](#run-the-docker-image-detached)
- [Password for superuser privileges](#password-for-superuser-privileges)

----------------------
#### Install Docker

Follow the instructions at the [Docker website](https://docs.docker.com/installation/)

Test Docker installation

    docker run hello-world

This command downloads a test image and runs it in a container.

----------------------
#### Create Docker image for grabb

There are two ways to do this:

* Pull the repository

        docker pull brankovics/grabb

* Create local docker repository

        git clone https://github.com/b-brankovics/grabb
        cd grabb/docker
        docker build -t localhost:5000/$USER/grabb .

Check if the image is created

    docker images

The image should be listed there

----------------------
#### Run the Docker image and attach local folders

To mount a folder to the Docker container specify the absolute
path on the host and the guest (inside the container). Using
the **-v** argument followed by the two paths
_host directory_**:**_guest directory_

Multiple folders maybe mounted to the same container.

If the **-rm** argument is added then the container will be
deleted after it is exited. (The Docker image won't be deleted)

**grabb** is the name of the Docker image that contains the desired installation of `GRAbB.pl`

    docker run -ti -v /absolute/path/to/the/directory/:/home/grabb/directory/ grabb

Now all the commands are issued are inside the container.
After the analysis is finished **Ctrl+D** has to be typed to return to the host computer.

----------------------
#### Run the Docker image detached

Same as for the previous example, it is possible to mount local folders.

Instead of using the **-ti** arguments, now the **-d** argument has to be used.

In addition, the `GRAbB.pl` command has to be specified as well. The output folder has
to be inside a directory that is mounted from the host computer,
because otherwise it is not possible to access the results.

    docker run -d -v /absolute/path/to/the/directory/:/home/grabb/directory/ grabb GRAbB.pl --ref <directory/ref.fas> --reads <directory/read*.fastq> --folder <derectory/folder> --prefix <prefix>

Docker will print the identifier of the container. Using this identifier it is possible to check on the run.

    docker logs <identifier>

At the end the results can be found in the folder specified as output in the host computer.

To remove the container:

    docker rm <identifier>

----------------------
#### Password for superuser privileges

The password is "**grabb**"