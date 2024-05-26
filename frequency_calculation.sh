#!/bin/bash
PROJECT_DIR="/home/projects/cge/data/projects/2024/projects_students/attila_beleon/thesis/data/hmmer_results"

awk '
    !/^#/ {  # Skip comment lines
        queryID = $4;  # Query ID
        targetID = $1;  # Target ID
        eValue = $7;  # E-value

        # Check if this target has not been seen or if this e-value is lower than previously seen for this target
        if (!seen[targetID, queryID] || eValue < evalue[targetID, queryID]) {
            evalue[targetID, queryID] = eValue;  # Store the lowest e-value for this target-query pair
            seen[targetID, queryID] = 1;  # Mark this target-query pair as seen
        }
    }
    END {
        # Count occurrences of the lowest e-value hits per target ID
        for (target_query in seen) {
            split(target_query, tq, SUBSEP);  # SUBSEP is the default subscript separator used by awk
            targetID = tq[1];
            queryID = tq[2];
            lowestHit[targetID]++;
        }
        # Print the target IDs and the number of times they were hit by the lowest e-value query
        for (targetID in lowestHit) {
            print targetID, lowestHit[targetID];
        }
    }
' "${PROJECT_DIR}/bla_AMRFinder_scores.domtblout" | sort -k2,2nr > "${PROJECT_DIR}/bla_AMRFinder_scores_freq.txt"
