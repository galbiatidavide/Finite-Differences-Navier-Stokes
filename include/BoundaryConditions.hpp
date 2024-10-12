#include <functional>
#include <string>

struct Component{
    Vec variable;
    DMStagStencilLocation location[2];
    std::string name;
};


class BoundaryConditions {

    private:
    PetscReal pi = 3.14159265358979323846;
    PetscReal eps = std::numeric_limits<float>::max();
    PetscReal A{4*sqrt(2)/(3*sqrt(3))};
    PetscReal vRef{4};
    PetscReal k{2*pi};
    PetscReal c{(5/6)*pi};
    PetscReal d{(1/6)*pi};
    PetscReal dt{0.01};
    
    protected:
    std::function<double(double, double, double, double)> bcFunctions[3];
    PetscReal iter;
    PetscReal theta = d*d*dt*iter;

    public: 

    BoundaryConditions() {
        bcFunctions[0] = [](double x, double y, double z, double theta) { return x + y + z;};
        bcFunctions[1] = [](double x, double y, double z, double theta) { return 2*x + 2*y + 2*z; };
        bcFunctions[2] = [](double x, double y, double z, double theta) { return 3*x+ 3*y + 3*z; };
    }

    PetscReal get_bcFunction(int index, double x, double y, double z, double theta) {
        if (index < 0 || index >= 3) {
            std::cerr << "Invalid boundary function index." << std::endl;
            return 0;
        }
        return bcFunctions[index](x, y, z, theta);
    }

    void set_bcFunction(int index, std::function<double(double, double, double, double)> func) {
        if (index < 0 || index >= 3) {
            std::cerr << "Invalid boundary function index." << std::endl;
            return;
        }
        bcFunctions[index] = func;
    }

    void set_time(){
        theta = d*d*dt*(iter+1);
    }

    PetscErrorCode set_IC(DM const & dmGrid, std::vector<Component> component){
    
        PetscFunctionBegin;

        for(unsigned int i = 0; i < component.size(); i++){

            PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3];
            DM              dmCoord;
            Vec             vecLocal, coord, coordLocal;
            PetscReal ****arrVec, ****arrCoord;

            DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
            DMGetCoordinateDM(dmGrid, &dmCoord);

            DMGetCoordinates(dmGrid, &coord);
            DMGetLocalVector(dmCoord, &coordLocal);
            DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

            for (int d = 0; d < 3; ++d) {
                DMStagGetLocationSlot(dmCoord, component[i].location[0], d, &icux[d]);
            }

            DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

            DMStagGetLocationSlot(dmGrid, component[i].location[0], 0, &iux);
            DMGetLocalVector(dmGrid, &vecLocal);
            DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

            for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
                for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
                    for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                        arrVec[ez][ey][ex][iux] = bcFunctions[i](arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
                    }
                }
            }

            DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
            DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
            DMLocalToGlobal(dmGrid, vecLocal, ADD_VALUES, component[i].variable);
            DMRestoreLocalVector(dmCoord, &coordLocal);
            DMRestoreLocalVector(dmGrid, &vecLocal);
        }

        PetscFunctionReturn(0);
    }

 PetscErrorCode set_BC(DM const & dmGrid, std::vector<Component> & components, std::array<PetscInt, 3> n_discr, Vec &globalVec){

    PetscFunctionBegin;
    
    for(unsigned int i = 0; i < components.size(); i++){

        PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iuxStart, iuxEnd, icuxStart[3], icuxEnd[3];
        DM              dmCoord;
        Vec             vecLocal, coord, coordLocal;
        PetscReal ****arrVec, ****arrCoord;

        DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (int d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, components[i].location[0], d, &icuxStart[d]);
            DMStagGetLocationSlot(dmCoord, components[i].location[1], d, &icuxEnd[d]);
        }

        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagGetLocationSlot(dmGrid, components[i].location[0], 0, &iuxStart);
        DMStagGetLocationSlot(dmGrid, components[i].location[1], 0, &iuxEnd);
        DMGetLocalVector(dmGrid, &vecLocal);
        DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

        // Set boundary conditions based on the component index
        if (i == 0) {
            // First component: Right and Left faces (ex == 0 for left, ex == start[0] + n[0] + nExtra[0] - 1 for right)
            for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
                for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                    if(start[0]==0)
                    arrVec[ez][ey][start[0]][iuxStart] = bcFunctions[i](arrCoord[ez][ey][start[0]][icuxStart[0]], arrCoord[ez][ey][start[0]][icuxStart[1]], arrCoord[ez][ey][start[0]][icuxStart[2]], theta);
                    if(start[0] + n[0] == n_discr[0])
                    arrVec[ez][ey][start[0] + n[0] - 1][iuxEnd] = bcFunctions[i](arrCoord[ez][ey][start[0] + n[0] - 1][icuxEnd[0]], arrCoord[ez][ey][start[0] + n[0] - 1][icuxEnd[1]], arrCoord[ez][ey][start[0] + n[0] - 1][icuxEnd[2]], theta);
                }
            }
        } else if (i == 1) {
            // Second component: Up and Down faces (ey == 0 for down, ey == start[1] + n[1] + nExtra[1] - 1 for up)
            for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
                for (ex = start[0]; ex < start[0] + n[0]; ++ex) {
                    if(start[1]==0)
                    arrVec[ez][start[1]][ex][iuxStart] = bcFunctions[i](arrCoord[ez][start[1]][ex][icuxStart[0]], arrCoord[ez][start[1]][ex][icuxStart[1]], arrCoord[ez][start[1]][ex][icuxStart[2]], theta);
                    if(start[1] + n[1] == n_discr[1])
                    arrVec[ez][start[1] + n[1] - 1][ex][iuxEnd] = bcFunctions[i](arrCoord[ez][start[1] + n[1] - 1][ex][icuxEnd[0]], arrCoord[ez][start[1] + n[1] - 1][ex][icuxEnd[1]], arrCoord[ez][start[1] + n[1] - 1][ex][icuxEnd[2]], theta);
                }
            }
        } else if (i == 2) {
            // Third component: Front and Back faces (ez == 0 for front, ez == start[2] + n[2] + nExtra[2] - 1 for back)
            for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
                for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                    if(start[2] == 0)
                    arrVec[start[2]][ey][ex][iuxStart] = bcFunctions[i](arrCoord[start[2]][ey][ex][icuxStart[0]], arrCoord[start[2]][ey][ex][icuxStart[1]], arrCoord[start[2]][ey][ex][icuxStart[2]], theta);
                    if(start[2] + n[2] == n_discr[2])
                    arrVec[start[2] + n[2] - 1][ey][ex][iuxEnd] = bcFunctions[i](arrCoord[start[2] + n[2] - 1][ey][ex][icuxEnd[0]], arrCoord[start[2] + n[2] - 1][ey][ex][icuxEnd[1]], arrCoord[start[2] + n[2] - 1][ey][ex][icuxEnd[2]], theta);
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
        DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, components[i].variable);
        DMRestoreLocalVector(dmCoord, &coordLocal);
        DMRestoreLocalVector(dmGrid, &vecLocal);
    }

    PetscFunctionReturn(0);
}


    ~BoundaryConditions() {};


};