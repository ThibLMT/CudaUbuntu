/**
 *\file allocate_variables.c
 *\brief body of function allocate_variables
 */

#include "deallocate_variables.h"
#include "def_global_variables.h"


void deallocate_variables(void)
{
  int ia,ja,ka;
/*
  // deallocation of array mic : It is not necessary because the memory is automatically released on program termination 
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  for (ka=0; ka<syssizez; ka++)
	    {
	      free(mic[ia][ja][ka]);
	      mic[ia][ja][ka]=NULL;
	    }
	  free(mic[ia][ja]);
	  mic[ia][ja]=NULL;
	}
      free(mic[ia]);
      mic[ia]=NULL;
    }     
  free(mic);
  mic=NULL;
	
  // deallocation of array mic_boundary :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  for (ka=0; ka<syssizez; ka++)
	    {
	      free(mic_boundary[ia][ja][ka]);
	      mic_boundary[ia][ja][ka]=NULL;
	    }
	  free(mic_boundary[ia][ja]);
	  mic_boundary[ia][ja]=NULL;
	}
      free(mic_boundary[ia]);
      mic_boundary[ia]=NULL;
    }     
  free(mic_boundary);
  mic_boundary=NULL;
	
  // deallocation of array mic_insert :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  free(mic_insert[ia][ja]);
	  mic_insert[ia][ja]=NULL;
	}
      free(mic_insert[ia]);
      mic_insert[ia]=NULL;
    }     
  free(mic_insert);
  mic_insert=NULL;
  
  // deallocation of array mic_insert_boundary :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  free(mic_insert_boundary[ia][ja]);
	  mic_insert_boundary[ia][ja]=NULL;
	}
      free(mic_insert_boundary[ia]);
      mic_insert_boundary[ia]=NULL;
    }     
  free(mic_insert_boundary);
  mic_insert_boundary=NULL;
  */
  // deallocation of the particle array
 free(particle);
	

} // end deallocate_variables


void deallocate_mic(void)
{
  int ia,ja,ka;

  // deallocation of array mic : It is not necessary because the memory is automatically released on program termination 
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  for (ka=0; ka<syssizez; ka++)
	    {
	      free(mic[ia][ja][ka]);
	      mic[ia][ja][ka]=NULL;
	    }
	  free(mic[ia][ja]);
	  mic[ia][ja]=NULL;
	}
      free(mic[ia]);
      mic[ia]=NULL;
    }     
  free(mic);
  mic=NULL;
	
  // deallocation of array mic_boundary :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  for (ka=0; ka<syssizez; ka++)
	    {
	      free(mic_boundary[ia][ja][ka]);
	      mic_boundary[ia][ja][ka]=NULL;
	    }
	  free(mic_boundary[ia][ja]);
	  mic_boundary[ia][ja]=NULL;
	}
      free(mic_boundary[ia]);
      mic_boundary[ia]=NULL;
    }     
  free(mic_boundary);
  mic_boundary=NULL;
	
  // deallocation of array mic_insert :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  free(mic_insert[ia][ja]);
	  mic_insert[ia][ja]=NULL;
	}
      free(mic_insert[ia]);
      mic_insert[ia]=NULL;
    }     
  free(mic_insert);
  mic_insert=NULL;
  
  // deallocation of array mic_insert_boundary :
  for (ia=0; ia<syssizex; ia++)
    {
      for (ja=0; ja<syssizey; ja++)
	{
	  free(mic_insert_boundary[ia][ja]);
	  mic_insert_boundary[ia][ja]=NULL;
	}
      free(mic_insert_boundary[ia]);
      mic_insert_boundary[ia]=NULL;
    }     
  free(mic_insert_boundary);
  mic_insert_boundary=NULL;
  
	

} // end deallocate_variables
