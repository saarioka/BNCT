source /home/ssaariok/software/geant4-v11.3.2-install/bin/geant4.sh
export G4PARTICLEHPDATA=/home/ssaariok/software/G4TENDL1.4

cmake ..
make -j20

ProtonLibraries=("G4TENDL1.4")
NeutronLibraries=("G4NDL4.7.1" "ENDF-VII.1" "ENDF-VIII.0" "BROND-3.1" "JEFF-3.3" "JENDL-4.0u")

for PLib in "${ProtonLibraries[@]}"; do
    for NLib in "${NeutronLibraries[@]}"; do
        export RUN_ID="${PLib}_${NLib}"
        export G4PROTONHPDATA=/home/ssaariok/software/${PLib}/Proton
        export G4NEUTRONHPDATA=/home/ssaariok/software/${NLib}
        echo -e "\nRunning with proton library: ${PLib}, neutron library: ${NLib}"
        echo "  RUN_ID: ${RUN_ID}"
        echo "  G4PROTONHPDATA: ${G4PROTONHPDATA}"
        echo "  G4NEUTRONHPDATA: ${G4NEUTRONHPDATA}"
        time ./exampleB2b run_proton.mac > "${RUN_ID}.log" 2>&1
        echo "Completed run for ${RUN_ID}"

        sleep 1
    done
done
