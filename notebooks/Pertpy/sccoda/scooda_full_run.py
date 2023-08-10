import os
from scooda_analysis import run_analysis

if __name__ == '__main__':

    for config in os.listdir("sccoda_analysis_input"):
        
        if not config.startswith("config"):
            continue

        print(f"Running analysis for {config}")
        run_analysis("sccoda_analysis_input/" + config)
        #run_analysis("sccoda_analysis_input/" + config, 0.4, 0.5)
        
    