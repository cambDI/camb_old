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

void R_standardiseMolecules(char **structures_file, char **standardised_file, int *isSDFInt, int *isTrainingInt, char **name_file, int *limit) {
    indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
    indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
   
    bool isSDF = (*isSDFInt!=0);
    bool isTraining = (*isTrainingInt!=0);
    bool useNamefile = (strcmp(*name_file,"")!=0);
    int l = *limit;
   
    bool debug = false;
   
    int structure, structureIter;
    int sdfWriter = indigoWriteFile(*standardised_file);
    ofstream namefile;
    if(useNamefile) namefile.open(*name_file);
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
    bool noEorW = true; // no errors or warnings
    while (structure = indigoNext(structureIter)) {
        readCount++;
        if(l != -1 && readCount > l) break;
        if (indigoCountAtoms(structure) != -1) {
            int structureIndex = indigoIndex(structure);
            
            //printf("Structure Index: %d, readCount: %d\n", structureIndex+1, readCount);
            string structureName = indigoName(structure);
            int structureClone = indigoClone(structure);
                
            if (debug) printf("folding hydrogens\n");
            indigoFoldHydrogens(structure);
            
            if((structureIndex+1)%50 == 0) printf(".");
            if((structureIndex+1)%5000 == 0) printf("\n");
            
            if (debug) printf("checking bad valence and ambiguousH\n");
            // skip over if indigo determines bad valence
            if( indigoCheckBadValence(structure)==NULL ) {
                Rprintf("%s (#%d) skipped over: INDIGO_BAD_VALANCE\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            // skip over if indigo determines ambiguous H
            if( indigoCheckAmbiguousH(structure)==NULL ) {
                Rprintf("%s (#%d) skipped over: INDIGO_AMBIGUOUSH\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            
            if (debug) printf("checking if organic\n");
            // skip over if not organic
            if(!isOrganic(structure)) {
                inorganicCount++;
                Rprintf("%s (#%d) skipped over: NOT_ORGANIC\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            
            bool wasAromatic = isAromatic(structure);
            // check if the structure can be dearomatized
            if (debug) printf("dearomitizing\n");
            indigoDearomatize(structure);
            
            if (debug) printf("checking aromatic\n");
            // if the structure is now not aromatic (i.e. hasn't failed the dearomitisation) then pass it through the inchi plugin
            if(!isAromatic(structure)) {
                if (debug) printf("passing through inchi\n");
                const char* inchi = indigoInchiGetInchi(structure);
                string trimmedInchi = trimInchi(inchi);
                int temp = indigoInchiLoadMolecule(trimmedInchi.c_str());
                if(strcmp(indigoInchiGetWarning(), "") != 0) {
                    Rprintf("%s (#%d) warning: while converting to Inchi: %s\n", structureName.c_str(), structureIndex+1, indigoInchiGetWarning());
                    noEorW = false;
                }
                else {
                    indigoFree(structure);
                    structure = temp;
                }
            }
            else {
                if(wasAromatic) {
                    Rprintf("%s (#%d) warning: failed dearomatize, didn't go through inchi\n", structureName.c_str(), structureIndex+1);
                    noEorW = false;
                }
            }
            
            // reset the isotopes of all the atoms in the molecule
            if (debug) printf("resetting isotopes\n");
            resetIsotopes(structure);
            
            // print warnings on molecular mass variations
            if (debug) printf("various conditions\n");
            float molecularMass = indigoMolecularWeight(structure);
            if(molecularMass<20) {
                tooLightCount++;
                Rprintf("%s (#%d) warning: molecular mass less than 20 daltons, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            if(molecularMass>900) {
                tooHeavyCount++;
                Rprintf("%s (#%d) warning: molecular mass greater than 900 daltons, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highFlourine = containsMoreThanX(structure, 3, 9);
            if(highFlourine) {
                highFlourineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 flourines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highChlorine = containsMoreThanX(structure, 3, 17);
            if(highChlorine) {
                highChlorineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 chlorines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highBromine = containsMoreThanX(structure, 3, 35);
            if(highBromine) {
                highBromineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 bromines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highIodine = containsMoreThanX(structure, 3, 53);
            if(highIodine) {
                highIodine++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 iodines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            
            // remove 'bad' molecules from the training set
            if(isTraining) {
                if(molecularMass < 20 || molecularMass > 900 || highFlourine || highChlorine || highBromine || highIodine) {
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue; 
                }
            }
            
            // modify the structure
            //if (debug) printf("neutralising\n");
            //int negCorrect = neutraliseCSPnegO(structure);
            //int posCorrect = neutraliseNposH(structure);
            
            // rendering example left in for reference
            //if(negCorrect > 0 || posCorrect > 0) {
            //    printf("structure index: %d\n", structureIndex+1);
            //    printf("negCorrect: %d PosCorrect: %d\n", negCorrect, posCorrect);
            //    renderPair(structureClone, structure, "render", "molecule", structureIndex+1, "Bad Molecule");
            //}
        
            // use the largest substructure
            if (debug) printf("selecting the largest substructure\n");
            int temp2 = pickLargestSubstructure(structure);
            indigoFree(structure);
            structure = temp2;
            
            // unfold the hydrogens (this doesn't happen anymore)
            //if (debug) printf("unfolding hydrogens\n");
            //indigoUnfoldHydrogens(structure);
            
            writeCount++;
            
            if(useNamefile) {
                namefile << std::string(indigoName(structureClone)) + "\n";
                char number[15];
                sprintf(number, "%d", writeCount);
                indigoSetName(structure, number);
            }
            else {
                indigoSetName(structure, indigoName(structureClone));
                Rprintf("%s\n", indigoName(structureClone));
            }
            indigoSdfAppend(sdfWriter, structure);
            indigoFree(structureClone);
            indigoFree(structure);
        }
        else {
            Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            noEorW = false;
        }
    }  
    printf("\n");
    indigoFree(structureIter);
    indigoClose(sdfWriter);
    indigoFree(sdfWriter);
    if(useNamefile) namefile.close();
   
    if(isTraining) {
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
   
    if(noEorW) {
        Rprintf("No errors or warnings encountered."); 
    }
}





















void R_standardiseMoleculesNoProtonation(char **structures_file, char **standardised_file, int *isSDFInt, int *isTrainingInt, char **name_file, int *limit) {
    indigoSetOption("aromaticity-model", "generic"); // enable the dearomatization to work properly
    indigoSetOption("dearomatize-verification", "false"); // enable the dearomatization to work properly
   
    bool isSDF = (*isSDFInt!=0);
    bool isTraining = (*isTrainingInt!=0);
    bool useNamefile = (strcmp(*name_file,"")!=0);
    int l = *limit;
   
    bool debug = false;
   
    int structure, structureIter;
    int sdfWriter = indigoWriteFile(*standardised_file);
    ofstream namefile;
    if(useNamefile) namefile.open(*name_file);
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
    bool noEorW = true; // no errors or warnings
    while (structure = indigoNext(structureIter)) {
        readCount++;
        if(l != -1 && readCount > l) break;
        if (indigoCountAtoms(structure) != -1) {
            int structureIndex = indigoIndex(structure);
            
            //printf("Structure Index: %d, readCount: %d\n", structureIndex+1, readCount);
            string structureName = indigoName(structure);
            int structureClone = indigoClone(structure);
                
            if (debug) printf("folding hydrogens\n");
            indigoFoldHydrogens(structure);
            
            if((structureIndex+1)%50 == 0) printf(".");
            if((structureIndex+1)%5000 == 0) printf("\n");
            
            if (debug) printf("checking bad valence and ambiguousH\n");
            // skip over if indigo determines bad valence
            if( indigoCheckBadValence(structure)==NULL ) {
                Rprintf("%s (#%d) skipped over: INDIGO_BAD_VALANCE\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            // skip over if indigo determines ambiguous H
            if( indigoCheckAmbiguousH(structure)==NULL ) {
                Rprintf("%s (#%d) skipped over: INDIGO_AMBIGUOUSH\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            
            if (debug) printf("checking if organic\n");
            // skip over if not organic
            if(!isOrganic(structure)) {
                inorganicCount++;
                Rprintf("%s (#%d) skipped over: NOT_ORGANIC\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
                indigoFree(structureClone);
                indigoFree(structure);
                continue;
            }
            
            bool wasAromatic = isAromatic(structure);
            // check if the structure can be dearomatized
            if (debug) printf("dearomitizing\n");
            indigoDearomatize(structure);
            
            if (debug) printf("checking aromatic\n");
            // if the structure is now not aromatic (i.e. hasn't failed the dearomitisation) then pass it through the inchi plugin
            if(!isAromatic(structure)) {
                if (debug) printf("passing through inchi\n");
                const char* inchi = indigoInchiGetInchi(structure);
                string trimmedInchi = trimInchi(inchi);
                int temp = indigoInchiLoadMolecule(trimmedInchi.c_str());
                if(strcmp(indigoInchiGetWarning(), "") != 0) {
                    Rprintf("%s (#%d) warning: while converting to Inchi: %s\n", structureName.c_str(), structureIndex+1, indigoInchiGetWarning());
                    noEorW = false;
                }
                else {
                    indigoFree(structure);
                    structure = temp;
                }
            }
            else {
                if(wasAromatic) {
                    Rprintf("%s (#%d) warning: failed dearomatize, didn't go through inchi\n", structureName.c_str(), structureIndex+1);
                    noEorW = false;
                }
            }
            
            // reset the isotopes of all the atoms in the molecule
            if (debug) printf("resetting isotopes\n");
            resetIsotopes(structure);
            
            // print warnings on molecular mass variations
            if (debug) printf("various conditions\n");
            float molecularMass = indigoMolecularWeight(structure);
            if(molecularMass<20) {
                tooLightCount++;
                Rprintf("%s (#%d) warning: molecular mass less than 20 daltons, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            if(molecularMass>900) {
                tooHeavyCount++;
                Rprintf("%s (#%d) warning: molecular mass greater than 900 daltons, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highFlourine = containsMoreThanX(structure, 3, 9);
            if(highFlourine) {
                highFlourineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 flourines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highChlorine = containsMoreThanX(structure, 3, 17);
            if(highChlorine) {
                highChlorineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 chlorines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highBromine = containsMoreThanX(structure, 3, 35);
            if(highBromine) {
                highBromineCount++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 bromines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            bool highIodine = containsMoreThanX(structure, 3, 53);
            if(highIodine) {
                highIodine++;
                Rprintf("%s (#%d) warning: molecule contains more than 3 iodines, may not give good prediction\n", structureName.c_str(), structureIndex+1);
                noEorW = false;
            }
            
            // remove 'bad' molecules from the training set
            if(isTraining) {
                if(molecularMass < 20 || molecularMass > 900 || highFlourine || highChlorine || highBromine || highIodine) {
                    indigoFree(structureClone);
                    indigoFree(structure);
                    continue; 
                }
            }
            
            // modify the structure
            //if (debug) printf("neutralising\n");
            //int negCorrect = neutraliseCSPnegO(structure);
            //int posCorrect = neutraliseNposH(structure);
            
            // rendering example left in for reference
            //if(negCorrect > 0 || posCorrect > 0) {
            //    printf("structure index: %d\n", structureIndex+1);
            //    printf("negCorrect: %d PosCorrect: %d\n", negCorrect, posCorrect);
            //    renderPair(structureClone, structure, "render", "molecule", structureIndex+1, "Bad Molecule");
            //}
        
            // use the largest substructure
            if (debug) printf("selecting the largest substructure\n");
            int temp2 = pickLargestSubstructure(structure);
            indigoFree(structure);
            structure = temp2;
            
            // unfold the hydrogens (this doesn't happen anymore)
            //if (debug) printf("unfolding hydrogens\n");
            //indigoUnfoldHydrogens(structure);
            
            writeCount++;
            
            
            if(useNamefile) {
                namefile << std::string(indigoName(structureClone)) + "\n";
                char number[15];
                sprintf(number, "%d", writeCount);
                indigoSetName(structure, number);
            }
            else {
                indigoSetName(structure, indigoName(structureClone));
                Rprintf("%s\n", indigoName(structureClone));
            }
            indigoSdfAppend(sdfWriter, structure);
            indigoFree(structureClone);
            indigoFree(structure);
        }
        else {
            Rprintf("%s, Readcount: %d\n", indigoGetLastError(), readCount);
            noEorW = false;
        }
    }  
    printf("\n");
    indigoFree(structureIter);
    indigoClose(sdfWriter);
    indigoFree(sdfWriter);
    if(useNamefile) namefile.close();
   
    if(isTraining) {
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
   
    if(noEorW) {
        Rprintf("No errors or warnings encountered."); 
    }
}

} // extern "C"







