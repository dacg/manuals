
from pylmgc90.chipy import *
import time


###### Funtions ####################
# ID for interactions: 1 (dkdk), 4 (dkpl), 5 (dkjc), 11 (plpl), 12(pljc)
def gap_values(id_int, tol_gap):
   # Verifying gap
   gap_max = 0
   gap_sum = 0
   gap_avg = 0
   counter_valid_inter = 0

   curr_handler = inter_handler_2D_getAll(id_int)

   #[0:1,icdan] = verlet_inter(icdtac)%coor(:,iadj)    # Contact point coordinates
   #[2:3,icdan] = verlet_inter(icdtac)%nuc(:,iadj)     # Local frame normal vector
   #[4  ,icdan] = verlet_inter(icdtac)%rln(iadj) / H   # Normal force
   #[5  ,icdan] = verlet_inter(icdtac)%rlt(iadj) / H   # Tangential force
   #[6  ,icdan] = verlet_inter(icdtac)%gaptt(iadj)     # Gap


   # Looping the contact list
   for i in range(0,len(curr_handler),1):
    # The current gap
    curr_gap = curr_handler[i,6]

    # Looking for the maximal gap (negative)
    if (gap_max > curr_gap): 
      gap_max=curr_gap
    # Computing the average gap
    if (curr_gap <= 0):
      counter_valid_inter += 1
      gap_sum += curr_gap
    #

    # If the average cannot be computed, lack of valid inters.
    if (counter_valid_inter>0):
      gap_avg = gap_sum/counter_valid_inter
    else :
      gap_avg = 0
    #

    #
    if (abs(gap_avg)>tol_gap):
      print ("Gap is too big!")
      CloseDisplayFiles()
      sys.exit()
    #
    #
   #print ("###############################")
   print ("###############################")
   print ("Average gap: " + str(gap_avg))
   print ("Maximal gap: " + str(gap_max))
   print ('valid_inter: '+ str(counter_valid_inter))
   #print ("###############################")
   #print ("###############################")
   return [gap_avg, gap_max, counter_valid_inter]
   
############################
Initialize()

checkDirectories()

# desactivation des messages de log
#utilities_DisableLogMes()

####
# info gestion du temps
dt = 5.e-4
theta = 0.5
nb_steps = 500000

# bavardage de certaines fonctions
echo = 0

# info generation fichier visu
freq_display = 1000. 
ref_radius   = 1000.

# parameter of check gap an criterial
frec_check     = 500
frec_check_gap = 1000
#############################################
d_mean = 4.99666888741

#############################################
tol_gap = 0.1*d_mean      ####
ids = [1,4,5,11,12]    # ID for interactions: 1 (dkdk), 4 (dkpl), 5 (dkjc), 11 (plpl), 12(pljc)

# info contact

#       123456789012345678901234567890
type = 'Stored_Delassus_Loops         '
norm = 'Quad '
tol = 0.1666e-3
relax = 1.0
gs_it1 = 81
gs_it2 = 201

SetDimension(2)
### definition des parametres du calcul ###
utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)

### lecture du modele ###

### model reading ###
utilities_logMes('READ BODIES')
RBDY2_ReadBodies()

utilities_logMes('READ INI DOF')
TimeEvolution_ReadIniDof()
RBDY2_ReadIniDof()

utilities_logMes('READ BEHAVIOURS')
bulk_behav_ReadBehaviours()
tact_behav_ReadBehaviours()

#LOADS
DISKx_LoadTactors()
POLYG_LoadTactors()
JONCx_LoadTactors()
RBDY2_LoadBehaviours()

utilities_logMes('READ INI Vloc Rloc')
TimeEvolution_ReadIniVlocRloc()
DKDKx_ReadIniVlocRloc()
DKJCx_ReadIniVlocRloc()
DKPLx_ReadIniVlocRloc() ###
PLPLx_ReadIniVlocRloc()
PLJCx_ReadIniVlocRloc()

utilities_logMes('READ DRIVEN DOF')
RBDY2_ReadDrivenDof()

### ecriture paranoiaque du modele ###
utilities_logMes('WRITE BODIES')
overall_WriteBodies()
RBDY2_WriteBodies()

utilities_logMes('WRITE BEHAVIOURS')
bulk_behav_WriteBehaviours()
tact_behav_WriteBehaviours()

utilities_logMes('WRITE DRIVEN DOF')
overall_WriteDrivenDof()
RBDY2_WriteDrivenDof()

#
Initialize()

### post2D ##
OpenDisplayFiles()

utilities_logMes('COMPUTE MASS')
RBDY2_ComputeMass()

nb_rbdy2 = RBDY2_GetNbRBDY2()

#######################################3
# The pressure on the box in Pa
pressure_box = 10000 #Pa

# Array to store the evolution of the volume of the box
criterion_stop = np.zeros(2)

#create a file with volumetric behavior 


for k in range(1, nb_steps + 1, 1):
   #
   utilities_logMes('itere : '+str(k))
   #
   utilities_logMes('INCREMENT STEP')
   TimeEvolution_IncrementStep()
   RBDY2_IncrementStep()

   utilities_logMes('DISPLAY TIMES')
   TimeEvolution_DisplayStep()

   utilities_logMes('COMPUTE Fext')
   RBDY2_ComputeFext()
   #ingresar una rampa de presion
   RBDY2_BiaxialLoading(nb_rbdy2-3,pressure_box,nb_rbdy2,pressure_box,nb_rbdy2-2,pressure_box,nb_rbdy2-1,pressure_box)

   utilities_logMes('COMPUTE Fint')
   RBDY2_ComputeBulk()

   utilities_logMes('COMPUTE Free Vlocy')
   RBDY2_ComputeFreeVelocity()
   #
   utilities_logMes('SELECT PROX TACTORS')
   overall_SelectProxTactors()
   DKJCx_SelectProxTactors()
   DKDKx_SelectProxTactors()
   DKPLx_SelectProxTactors() #####
   PLPLx_SelectProxTactors()
   PLJCx_SelectProxTactors()

   #
   DKJCx_RecupRloc()
   DKDKx_RecupRloc()
   DKPLx_RecupRloc() ####
   PLPLx_RecupRloc()
   PLJCx_RecupRloc()
   nlgs_ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
   DKJCx_StockRloc()
   DKDKx_StockRloc()
   DKPLx_StockRloc() ####
   PLPLx_StockRloc()
   PLJCx_StockRloc()
   #
   utilities_logMes('COMPUTE DOF')
   RBDY2_ComputeDof()
   #
   utilities_logMes('UPDATE DOF')
   TimeEvolution_UpdateStep()
   RBDY2_UpdateDof()
   #
   utilities_logMes('WRITE LAST DOF')
   TimeEvolution_WriteLastDof()
   RBDY2_WriteLastDof()
   #
   utilities_logMes('WRITE LAST Vloc Rloc')
   TimeEvolution_WriteLastVlocRloc()
   DKDKx_WriteLastVlocRloc()
   DKJCx_WriteLastVlocRloc()
   DKPLx_WriteLastVlocRloc() ####
   PLPLx_WriteLastVlocRloc()
   PLJCx_WriteLastVlocRloc()
   #
   ### post2D ###
   WriteDisplayFiles(freq_display,ref_radius)

   ### wrtieout handling ###
   overall_CleanWriteOutFlags()

   # Adding a criterion to stop the simulation
   # This is a criterion in Volume. The volume should not change over a given 
   # tolerance in order to stop the iterations.

   # The reference
   if (k==1):
      Coor_d=RBDY2_GetBodyVector('Coor_',nb_rbdy2-3)
      Coor_u=RBDY2_GetBodyVector('Coor_',nb_rbdy2-2)
      Coor_l=RBDY2_GetBodyVector('Coor_',nb_rbdy2-1)
      Coor_r=RBDY2_GetBodyVector('Coor_',nb_rbdy2)

      # Computing the size of the box
      heigh_0 = abs(Coor_d[1] - Coor_u[1])
      length_0 = abs(Coor_l[0] - Coor_r[0])

      # Storing the volume of the box
      criterion_stop[0] = heigh_0*length_0

   # The normal cases
   else :
    if (k%frec_check == 0 and k >= 10000):  #revision
      Coor_d=RBDY2_GetBodyVector('Coor_',nb_rbdy2-3)
      Coor_u=RBDY2_GetBodyVector('Coor_',nb_rbdy2-2)
      Coor_l=RBDY2_GetBodyVector('Coor_',nb_rbdy2-1)
      Coor_r=RBDY2_GetBodyVector('Coor_',nb_rbdy2)

      heigh = abs(Coor_d[1] - Coor_u[1])
      length = abs(Coor_l[0] - Coor_r[0])

      # Storing the current volume
      criterion_stop[1] = heigh*length

      # I suppose the volume is always being reduced. I.E. T0-T1 always positive
      # If this is not the case, the box should be rebounding... Weird. 
      # So adding a warning and key entry to verify.

      if ((criterion_stop[0] - criterion_stop[1]) < 0.):
         print ("It's time to reset speeds to zero!!")
         break
      else : 
         def_vol = criterion_stop[1] - heigh_0*length_0
         f = open('file.txt', 'a')
         f.write(str(def_vol)+','+str(k*dt)+'\n')
         f.close()
         criterion_stop[0] = criterion_stop[1]

      #
    if  k%frec_check_gap == 0:  #revision
      # verifying gap 
      gap_max_t = 0
      gap_avg_t = 0
      gap_sum   = 0
      cont      = 0
      for j in range(0,len(ids),1):
        [gap_avg, gap_max, valid_inter] = gap_values(ids[j], tol_gap)
        
        if (gap_max_t > gap_max): 
          gap_max_t = gap_max
        
        if valid_inter > 0 :
          gap_sum += gap_avg
          cont += 1
        #

      # If the average cannot be computed, lack of valid inters.
      if (gap_sum<0):
        gap_avg_t = gap_sum/cont
      else :
        gap_avg_t = 0
    #
      print ("###############################")
      print ("Average gap: " + str(gap_avg_t))
      print ("Maximal gap: " + str(gap_max_t))
    #
      if (abs(gap_avg)>tol_gap):
         print ("Gap is too big!")
         CloseDisplayFiles()
         sys.exit()

      
      #
   #
#

CloseDisplayFiles()

Finalize()
