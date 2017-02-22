# YinMem
This is the repository for YinMem project: a parallel distributed indexed in-memory computing system

Components for YinMem:
  a. Hadoop cluster 
  b. Zookeeper 
  c. Accumulo database 
  d. Alluxio: in-memory file system 
  e. pMatlab. 
  
  In our experiment, we have set up a 32-nodes Hadoop cluster, n101-n132, n117 as the namenode. 
  
  Zookeeper nodes: n117,n118,n119.
  
  Accumulo master node: n117, others are workers.
  
  For alluxio, we deploy it in a local mode: 
  export ALLUXIO_UNDERFS_ADDRESS=${ALLUXIO_UNDERFS_ADDRESS:-${ALLUXIO_HOME}/underFSStorage}
  That means the underFS storage will be based on local HDD. 
  Also n117 is the master node 
  
  
  For pMatlab, more info can be found in the following link:
  http://www.ll.mit.edu/mission/cybersec/softwaretools/pmatlab/pmatlab.html
  
  
