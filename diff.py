def find_the_drift():
    # Adjust filenames if needed
    f1_path = 'blist7.txt'
    f2_path = 'johansson_stripped.txt'
    
    with open(f1_path, 'r') as f1, open(f2_path, 'r') as f2:
        # We only care about the length of the smaller file
        index = 0
        while True:
            char1 = f1.read(1)
            char2 = f2.read(1)
            
            if not char1: # End of 1M file
                print("\n✅ REACHED THE END: generated file is 100% matched to compareto file")
                return
            
            if char1 != char2:
                print(f"\n❌ DRIFT DETECTED at digit index: {index:,}")
                
                # Rewind slightly to show context
                f1.seek(max(0, index - 20))
                f2.seek(max(0, index - 20))
                
                context1 = f1.read(50)
                context2 = f2.read(50)
                
                print(f"Context (Index {max(0, index-20)} to {max(0, index+30)}):")
                print(f"generated File: ...{context1}...")
                print(f"compareto File: ...{context2}...")
                print(f"            {' ' * 20}^ DRIFT STARTS HERE")
                return
            
            index += 1
            if index % 2500 == 0:
                print(f"Verified {index:,} digits...")

find_the_drift()