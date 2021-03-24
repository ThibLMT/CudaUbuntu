//
// Created by ThibLMT on 23/03/2021.
//

#include "init_params.cuh"

void adi_params(material_data *prop_mat_part,geom_struct *geom)
{
    material_data tab_mat;
    tab_mat=*prop_mat_part;
    // E Young modulus Pa (N.m-2 --> N.unity-2)
    tab_mat.E=tab_mat.E*geom->unity*geom->unity;
    // Bulk density
    tab_mat.density=tab_mat.density*geom->unity*geom->unity*geom->unity;
    *prop_mat_part=tab_mat;
}

int get_line_char(char *detect,char **test)
{
    char *Nparafile;
    FILE *fparafile;
    static const char separateur[] = " ,()=\t\n[]";
    int nbjalon = 0;
    char *s;
    char ligne[256];
    char *jalon[256];
    int i = 0;
    int nstring,lstring;

    Nparafile="params.in";

    fparafile = fopen (Nparafile, "r");
    if(fparafile != NULL)
    {
        while(!feof (fparafile))
        {
            fgets (ligne, 256, fparafile);
            if (feof (fparafile)) break;
            // First call of strtok function
            s = (char *) strtok (ligne, (const char *) separateur);
            nbjalon = 0;
            jalon[0] = s;
            while (s != NULL)
            {
                jalon[nbjalon++] = s;
                s = (char *) strtok (NULL, (const char *) separateur);
            }
            jalon[nbjalon] = 0;

            // analyse de la ligne
            i = 0;
            if (!strcmp ((const char *) jalon[i],detect))
            {
                while (i < nbjalon)
                {
                    if ((!strncmp ((const char *) jalon[i], "#", 1))||(!strncmp ((const char *) jalon[i], "!", 1))||(!strncmp ((const char *) jalon[i], ";", 1)))
                    {nstring=i;break;}
                    i++;
                }
                break;
            }


        }

        fclose (fparafile);
    }
    else
    {
        fprintf(stderr,"Error! No parameter file (%s) \n",Nparafile);
        exit(errno);
    }
    for(i=0;i<nstring;i++){strcpy(test[i],jalon[i]);}
    return nstring + 1;
}

void read_table_mat(material_data *tab_mat_part)
{
    material_data mat_params;
    char **tab_char;  //!< array of string (50 elements)
    int i,k,nstring;

    // Initialization and allocation
    mat_params=*tab_mat_part;
    tab_char= static_cast<char **>(malloc(50 * sizeof(char *)));
    for(i=0;i<50;i++) {tab_char[i]= static_cast<char *>(malloc(256 * sizeof(char)));}
    // Read string "tab_mat" and grep info
    nstring=get_line_char("tab_mat",tab_char);

    for(k=0;k<nstring;k++)
    {printf(" %d  %s \n",k,tab_char[k]);
        if (!strcmp ((const char *) tab_char[k], "density"))
            mat_params.density = atof (tab_char[k + 1]);
        if (!strcmp ((const char *) tab_char[k], "E"))
            mat_params.E = atof (tab_char[k + 1]);
        if (!strcmp ((const char *) tab_char[k], "nu"))
            mat_params.nu = atof (tab_char[k + 1]);
        if (!strcmp ((const char *) tab_char[k], "cn"))
            mat_params.cn = atof (tab_char[k + 1]);
    }

    *tab_mat_part=mat_params;
}