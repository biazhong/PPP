# PPP
On this site, we include all the codes used to conduct the numerical experiments of the paper entitled "Speeding Up Paulson's Procedure for Large-Scale Problem Using Parallel Computing". In Section 7.1 of the paper, we conduct the numerical experiments in a single processor computing environment. All the codes are written in Java. In Section 7.2 of the paper, we conduct the numerical experiments in different parallel computing environments. All the codes are written in C++. To compile and run the Java codes, one may download integrated development environment (IDE) Eclipse from https://www.eclipse.org/. To compile and run the C++ codes, one may download parallel computing platform MPICH-3.3.2 we used in the paper from https://www.mpich.org/. Site https://mpitutorial.com/tutorials/ includes the tutorials on how to set up the parallel computing environments for MPI. For the numerical experiments conducted in Section 7.1, all the samples are drawn from gaussian distributions. Therefore, java.util.Random.nextGaussian() method is used to generated observations for different alternatives. For the numerical experiments conducted in Section 7.2, the simulation optimization problem we considered is the three-stage buffer allocation problem. The detailed description of the problem can be found at http://simopt.org. The functions used to generate observations for different alternatives are embeded in each procedure. All the C++ files in directory Section 7.2 are self-contained and  can be compiled and run independently. 


---Directory "Code for Section 7.1" contains all the Java codes that are used to conduct the numerical experiments listed in Section 7.1 of the paper entitled "Speeding Up Paulson's Procedure for Large-Scale Problems Using Parallel Computing". 

Files KN.java, ModifiedPaulson.java, and Paulson.java are the classes that implement the KN procedure, the PPP procedure, and the Paulson's procedure respectively in Section 7.1.

Files mainKN.java, mainModifiedPaulson.java, and mainPaulson.java are the main classes that run the experiments related to the KN procedure, the PPP proceudre, and the Paulson's procedure respectively in Section 7.1.

---Directory "Code for Section 7.2" contains all the codes that are used to conduct the numerical experiments listed in Section 7.2 of the paper entitled "Speeding Up Paulson's Procedure for Large-Scale Problems Using Parallel Computing".
  
  File EtaFunc.java is the Java codes that can be used to calculate the parameter eta used in the GSP procedure;
  
  File Rinott.java is the Java codes that can be used to calculate the Rinott's constant used in the GSP procedure;
  
  File GSP.cpp is the C++ codes that are used to conduct the experiments related to the GSP procedure in Section 7.2;
  
  File PPP-PACNew.cpp is the C++ codes that are used to conducted the experiments related to the PAC-PPP procedure in Tables 3, 5, 6 and Figure 6;
  
  File PPP-PACNewRandomSleep.cpp is the C++ codes that are used to conducted the experiment related to the PAC-PPP procedure in Table 4;
  
  File PPPNew.cpp is the C++ codes that are used to conducted the experiment related to the PPP procedure in Table 3;
  
  File PSS.cpp is the C++ codes that are used to conducted the experiment related to the PSS procedure in Table 3;
  
  File VKN.cpp is the C++ code that are used to conducted the experiment related to the VKN procedure in Table 3;
  
  File \_batchingc++RandomSleepTime.cpp contains the C++ codes used to conducted the experiment related to the batching method in Table 2;
  
  File \_onebyonec++RandomSleepTime.cpp contains the C++ codes used to conducted the experiment related to the one-by-one method in Table 2;


---The csv file "Numerical Results.csv" contains all the numerical results for the experiments we have conducted for the paper entitled "Speeding Up Paulson's Procedure for Large-Scale Problems Using Parallel Computing". Besides the results listed in the form of tables, it includes the exact numerical results for the experiments whose results are displayed in figures.



#NOTE: In directory "Code for Section 7.2", while implementing different procedures, we highlight the parameters that should be changed from one problem instance to another by commenting the parameter with "Input Parameter:..."

#NOTE: The exact computing environments used to conduct different experiments are clearly specified in the paper.

#NOTE: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability,
clarity, generality, and efficiency has not been well considered.
