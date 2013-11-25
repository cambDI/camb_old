#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "indigo.h"
#include "indigo-renderer.h"
#include "indigo-inchi.h"
#include "./smindigo.h"

using namespace std;

string trimInchi(const char* inchi) {
    string trimmed(inchi);
    size_t found;
    
    found = trimmed.find("/t");
    if(found != string::npos) {
        trimmed = trimmed.substr(0, found);
    }
    
    found = trimmed.find("/b");
    if(found != string::npos) {
        trimmed = trimmed.substr(0, found);
    }
    return trimmed;
}

extern "C" {
    
    void R_drawMoleculeInSDF(char **structures_file, int *structureNumberIn, char **filename, int *useNameAsTitleInt) {
        int structure, fileIter;
        int structureNumber = structureNumberIn[0];
        printf("structureNumber: %d\n", structureNumber);
        bool useNameAsTitle = (*useNameAsTitleInt!=0);
        fileIter = indigoIterateSDFile(*structures_file);
        int i = 0;
        while ((structure = indigoNext(fileIter))) {
            i++;
            if(i==structureNumber) {
                indigoFoldHydrogens(structure);
                printf("filename: %s", *filename);
                if(useNameAsTitle) {
                    renderMolecule(structure, *filename);
                }
                else {
                    renderMolecule(structure, "", *filename, -1, "");
                }
                break;
            }
            indigoFree(structure);
        }
        indigoFree(fileIter);
    }
    
    void R_drawMoleculeInSDFbyID(char **structures_file, char **structureIDIn, char **filename, int *useNameAsTitleInt) {
        int structure, fileIter;
        bool useNameAsTitle = (*useNameAsTitleInt!=0);
        fileIter = indigoIterateSDFile(*structures_file);
        int i = 0;
        while ((structure = indigoNext(fileIter))) {
            i++;
            char str[256];
            sprintf(str, "%s", *structureIDIn);
            if(strcmp(str,indigoName(structure)) == 0) {
                indigoFoldHydrogens(structure);
                if(useNameAsTitle) {
                    renderMolecule(structure, *filename);
                }
                else {
                    renderMolecule(structure, "", *filename, -1, "");
                }
                break;
            }
            indigoFree(structure);
        }
        indigoFree(fileIter);
    }
    
    void R_standardiseMolecules(char **structures_file,
                                char **standardised_file,
                                char **removed_file,
                               // char **target_field_name,
                                int *isSDFInt,
                                int *removeInorganicInt,
                                int *fluorineLimitInt,
                                int *chlorineLimitInt,
                                int *bromineLimitInt,
                                int *iodineLimitInt,
                                int *minMassLimitInt,
                                int *maxMassLimitInt,
                                int *numberProcessedInt) {
        indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
        indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
        
        bool isSDF = (*isSDFInt!=0);
        //bool write_targets = !(*target_field_name && *target_field_name[0] == '\0');
        bool removeInorganic = (*removeInorganicInt!=0);
        int numberProcessed = *numberProcessedInt;
        int fluorineLimit = *fluorineLimitInt;
        int chlorineLimit = *chlorineLimitInt;
        int bromineLimit = *bromineLimitInt;
        int iodineLimit = *iodineLimitInt;
        int minMassLimit = *minMassLimitInt;
        int maxMassLimit = *maxMassLimitInt;
        
        bool debug = false;
        
        int structure, structureIter;
        int sdfWriter = indigoWriteFile(*standardised_file);
        int removedWriter = indigoWriteFile(*removed_file);
        
        ofstream target_stream;
        
        //if(write_targets) {
            target_stream.open("targets.csv", ios::out | ios::trunc); // delete the current file
            target_stream.close();
            target_stream.open("targets.csv", ios::out | ios::ate | ios::app | ios::binary);
       // }
        
        if(isSDF) {
            Rprintf("Reading SDF (C)\n");
            structureIter = indigoIterateSDFile(*structures_file);
        }
        else {
            Rprintf("Reading SMILES (C)\n");
            structureIter = indigoIterateSmilesFile(*structures_file);
        }
        
        
        
        int readCount = 0;
        int writeCount = 0;
        int inorganicCount = 0;
        int tooLightCount = 0;
        int tooHeavyCount = 0;
        int highFlourineCount = 0;
        int highChlorineCount = 0;
        int highBromineCount = 0;
        int highIodineCount = 0;
        
        int props,prop,propFirst,propsFirst;
        int propss = 0;
        int namesprop = 0;
        while (structure = indigoNext(structureIter)) {
            readCount++;
            if(numberProcessed != -1 && readCount > numberProcessed) break;
            if (indigoCountAtoms(structure) != -1) {
                int structureIndex = indigoIndex(structure);
                double target_value = 0;
               // if(write_targets) {
                    
                    if(namesprop == 0){
                        propsFirst = indigoIterateProperties(structure);
                        while (propFirst = indigoNext(propsFirst)) {
                            string prop_val_first = indigoName(propFirst);
                            target_stream << prop_val_first << ",";
                            indigoFree(propFirst);
                        }
                        target_stream << "Kept" << "\n";
                        namesprop=1;
                    }
                    
                  //  if(indigoHasProperty(structure, *target_field_name)) {
                        ///target_value = atof(indigoGetProperty(structure, *target_field_name));
                        
                        props = indigoIterateProperties(structure);
                        while (prop = indigoNext(props)) {
                            
                            string target_value = indigoGetProperty(structure, indigoName(prop));
                            
                            target_stream << target_value << ",";
                            
                            indigoFree(prop);
                        }

 
                   // }
                    //else {
                    ///    write_targets = false;
                     ///   target_stream << "no property found on first molecule in file with fieldname: " << *target_field_name << endl;
                     ///   target_stream << "The molecules will not be standardised" << endl; /////modified

                        
                     ///   break; //// modified
                    //}
                //   }
                
                //printf("Structure Index: %d, readCount: %d\n", structureIndex+1, readCount);
                string structureName = indigoName(structure);
                int structureClone = indigoClone(structure);
                
                if (debug) Rprintf("%s (#%d)\n", structureName.c_str(), structureIndex+1);
                if (debug) Rprintf("folding hydrogens\n");
                indigoFoldHydrogens(structure);
                
                if((structureIndex+1)%50 == 0) printf(".");
                if((structureIndex+1)%5000 == 0) printf("\n");
                
                if (debug) Rprintf("checking bad valence and ambiguousH\n");
                // skip over if indigo determines bad valence
                if( indigoCheckBadValence(structure)==NULL ) {
                    Rprintf("%s (#%d) skipped over: INDIGO_BAD_VALANCE\n", structureName.c_str(), structureIndex+1);
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue;
                }
                // skip over if indigo determines ambiguous H
                if( indigoCheckAmbiguousH(structure)==NULL ) {
                    Rprintf("%s (#%d) skipped over: INDIGO_AMBIGUOUSH\n", structureName.c_str(), structureIndex+1);
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue;
                }
                
                if (debug) Rprintf("checking if organic\n");
                // skip over if not organic
                if(!isOrganic(structure)) {
                    inorganicCount++;
                    if(removeInorganic) {
                        Rprintf("%s (#%d) skipped over: NOT_ORGANIC\n", structureName.c_str(), structureIndex+1);
                        indigoFree(structureClone);
                        indigoFree(structure);
                        continue;
                    }
                }
                
                bool wasAromatic = isAromatic(structure);
                // check if the structure can be dearomatized
                if (debug) Rprintf("dearomitizing\n");
                indigoDearomatize(structure);
                
                if (debug) Rprintf("checking aromatic\n");
                // if the structure is now not aromatic (i.e. hasn't failed the dearomitisation) then pass it through the inchi plugin
                if(!isAromatic(structure)) {
                    if (debug) Rprintf("passing through inchi\n");
                    const char* inchi = indigoInchiGetInchi(structure);
                    string trimmedInchi = trimInchi(inchi);
                    int temp = indigoInchiLoadMolecule(trimmedInchi.c_str());
                    if(strcmp(indigoInchiGetWarning(), "") != 0) {
                        Rprintf("%s (#%d) warning: while converting to Inchi: %s\n", structureName.c_str(), structureIndex+1, indigoInchiGetWarning());
                    }
                    else {
                        indigoFree(structure);
                        structure = temp;
                    }
                }
                else {
                    if(wasAromatic) {
                        Rprintf("%s (#%d) warning: failed dearomatize, didn't go through inchi\n", structureName.c_str(), structureIndex+1);
                    }
                }

                
                // reset the isotopes of all the atoms in the molecule
                if (debug) Rprintf("resetting isotopes\n");
                resetIsotopes(structure);
                
                // print warnings on molecular mass variations
                if (debug) Rprintf("various conditions\n");
                float molecularMass = indigoMolecularWeight(structure);
                if(minMassLimit != -1 && molecularMass < minMassLimit) {
                    tooLightCount++;
                    Rprintf("%s (#%d) warning: molecular mass less than %d daltons\n", structureName.c_str(), structureIndex+1, minMassLimit);
                }
                if(maxMassLimit != -1 && molecularMass > maxMassLimit) {
                    tooHeavyCount++;
                    Rprintf("%s (#%d) warning: molecular mass greater than %d daltons\n", structureName.c_str(), structureIndex+1, maxMassLimit);
                }
                bool highFlourine = fluorineLimit != -1 && containsMoreThanX(structure, fluorineLimit, 9);
                if(highFlourine) {
                    highFlourineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d flourines\n", structureName.c_str(), structureIndex+1, fluorineLimit);
                }
                bool highChlorine = chlorineLimit != -1 && containsMoreThanX(structure, chlorineLimit, 17);
                if(highChlorine) {
                    highChlorineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d chlorines\n", structureName.c_str(), structureIndex+1, chlorineLimit);
                }
                bool highBromine = bromineLimit != -1 && containsMoreThanX(structure, bromineLimit, 35);
                if(highBromine) {
                    highBromineCount++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d bromines\n", structureName.c_str(), structureIndex+1, bromineLimit);
                }
                bool highIodine = iodineLimit != -1 && containsMoreThanX(structure, iodineLimit, 53);
                if(highIodine) {
                    highIodine++;
                    Rprintf("%s (#%d) warning: molecule contains more than %d iodines\n", structureName.c_str(), structureIndex+1, iodineLimit);
                }
                
                // remove 'bad' molecules (normally used during training)
                if(molecularMass < minMassLimit || molecularMass > maxMassLimit || highFlourine || highChlorine || highBromine || highIodine) {
                    indigoSdfAppend(removedWriter, structure);
                    indigoFree(structureClone);
                    indigoFree(structure);
                    target_stream << "0" << endl;
                    continue;
                } else {

                        target_stream << "1" << endl;

                }
                
                
                // use the largest substructure
                if (debug) Rprintf("selecting the largest substructure\n");
                int temp2 = pickLargestSubstructure(structure);
                indigoFree(structure);
                structure = temp2;
                
                writeCount++;
                std::string addition = "Standardised_";
                std::string name = indigoName(structureClone);
                indigoSetName(structure, (addition + name).c_str());
                //indigoSetProperty(structure, (addition + name).c_str(), (addition + name).c_str());
                indigoSdfAppend(sdfWriter, structure);
                

                indigoFree(structureClone);
                indigoFree(structure);
            }
            else {
                Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            }
        }  
        printf("\n");
        indigoFree(structureIter);
        indigoClose(sdfWriter);
        indigoFree(sdfWriter);
        indigoClose(removedWriter);
        indigoFree(removedWriter);
        target_stream.close();
        
        printf("\ninorganicCount: %d\n", inorganicCount);
        printf("tooLightCount: %d\n", tooLightCount);
        printf("tooHeavyCount: %d\n", tooHeavyCount);
        printf("highFlourineCount: %d\n", highFlourineCount);
        printf("highChlorineCount: %d\n", highChlorineCount);
        printf("highBromineCount: %d\n", highBromineCount);
        printf("highIodineCount: %d\n", highIodineCount);
        printf("readCount: %d\n", readCount);
        printf("writeCount: %d\n", writeCount);
    }
    
} // extern "C"







