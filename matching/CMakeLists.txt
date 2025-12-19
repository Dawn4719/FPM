set(SUBGRAPH_MATCHING_FILE_LIST
        matchingcommand.cpp
        matchingcommand.h
        FilterVertices.cpp
        FilterVertices.h
        BuildTable.cpp
        BuildTable.h
        GenerateQueryPlan.cpp
        GenerateQueryPlan.h
        EvaluateQuery.cpp
        EvaluateQuery.h
        GenerateFilteringPlan.cpp
        GenerateFilteringPlan.h
        FPM.cpp)

set(SUBGRAPH_MATCHING_FILE_LIST2
        matchingcommand.cpp
        matchingcommand.h
        FilterVertices.cpp
        FilterVertices.h
        BuildTable.cpp
        BuildTable.h
        GenerateQueryPlan.cpp
        GenerateQueryPlan.h
        EvaluateQuery.cpp
        EvaluateQuery.h
        GenerateFilteringPlan.cpp
        GenerateFilteringPlan.h
        FPMBase.cpp)

set(GRAPH_CONVERTER_FILE_LIST
        GraphConverter.cpp)

add_executable(FPM.out ${SUBGRAPH_MATCHING_FILE_LIST})

target_link_libraries(FPM.out PUBLIC graph utility )

add_executable(FPMbase.out ${SUBGRAPH_MATCHING_FILE_LIST2})

target_link_libraries(FPMbase.out PUBLIC graph utility )

#add_executable(SubgraphMatchingP.out ${SUBGRAPH_MATCHING_FILE_LIST3})

#target_link_libraries(SubgraphMatchingP.out PUBLIC graph utility )

#add_executable(GraphConverter.out ${GRAPH_CONVERTER_FILE_LIST})
#target_link_libraries(GraphConverter.out graph utility)
#
#add_executable(EdgeListConverter.out EdgeListToCSR.cpp)
#add_executable(ReassignLabel ReassignLabel.cpp)
#add_executable(SelectEdge SelectEdge.cpp)