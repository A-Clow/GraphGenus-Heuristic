from sage.all import *
import copy
import random
import itertools as it
import math


def random_embedding_for_sage(G):
    """
    Input: A sagemath graph object G
    Output: A randomly selected 2-cell embedding of G, recorded as an orientation system
    """
    
    embedding = {}
    
    for u in G.vertices():
        local_orientation = G.neighbors(u)
        random.shuffle(local_orientation)
        
        embedding[u] = local_orientation
        
    return embedding


def check_genus(G,fixed_embedding):
    """
    Input: A sagemath graph object G, and an embedding of G
    Output: The non-oriented genus of the embedding
    """

    n = G.order()
    e = G.size()
    f = len(list(G.faces(embedding=fixed_embedding)))

    return -((n-e+f-2)/2)


def search_for_low_genus(G,generations, k,current_generation_embbeding):
    
    """
    Input: A sagemath graph object G, and an embedding of G
    Output: Runs a genetic algorithm searching for a low genus embedding, optimising over k randomly sampled "columns", 
           with objective function the number of polyhedra faces, the search runs n iterations/generations. Returns the best embedding found.
    """

############# Be careful, this function runs in more than (maxDegree!)^k time (when max degree is small and k is small this is fine)
    
    best_embedding = current_generation_embbeding 
    counter = 0
    
    best_genus = check_genus(G,best_embedding)

    while counter < generations:

        #if stopper == True:
            #break

        counter = counter + 1
        Samples = []
        
        for i in range(10):
        
            sample  = random.sample(G.vertices(), k)
            Samples.append(sample)

        for vertex_sample in Samples:
            
            optimisation_dict = {}
            vector_list = []
            embeddings_to_check = []

            for u in vertex_sample:
                Options = list(it.permutations(G.neighbors(u)))
                optimisation_dict[u] = list(Options)
                vector_list.append(list(range(math.factorial(len(G.neighbors(u))))))

            vectors = list(it.product(*vector_list))
                

            for vec in vectors:
                
                new_embedding = copy.copy(best_embedding)
                index = 0

                while index < k:

                    new_embedding[vertex_sample[index]] = list(optimisation_dict[vertex_sample[index]][vec[index]])
                    index = index + 1
                    

                embeddings_to_check.append(new_embedding)

        
        for embedding in embeddings_to_check:
            new_genus = check_genus(G,embedding)

            if new_genus < best_genus:
                best_embedding  = embedding
                best_genus = new_genus
                

    print(best_genus)
    return best_embedding



def REGsearch_for_low_genus(G,generations, k,current_generation_embbeding):
    
    """
    Input: A sagemath regular graph object G, and an embedding of G
    Output: Runs a genetic algorithm searching for a low genus embedding, optimising over k randomly sampled "columns", 
           with objective function the number of polyhedra faces, the search runs n iterations/generations. Returns the best embedding found.
    """

############# Slightly optimized for regular graphs
    
    best_embedding = current_generation_embbeding 
    counter = 0

    # There is a bug latter if the graph is not regular, needs to be fixed later
    D = max(G.degree())
    d = min(G.degree())

    if d!= D:
        print("Error: please give me a regular graph")
        return
        
    vectors = list(it.product(range(math.factorial(d)), repeat = k))
    #
    
    best_genus = check_genus(G,best_embedding)

    while counter < generations:

        #if stopper == True:
            #break

        counter = counter + 1
        Samples = []
        
        for i in range(10):
        
            sample  = random.sample(G.vertices(), k)
            Samples.append(sample)

        for vertex_sample in Samples:
            
            optimisation_dict = {}
            embeddings_to_check = []

            for u in vertex_sample:
                Options = list(it.permutations(G.neighbors(u)))
                optimisation_dict[u] = list(Options)
                

            for vec in vectors:
                
                new_embedding = copy.copy(best_embedding)
                index = 0

                while index < k:

                    new_embedding[vertex_sample[index]] = list(optimisation_dict[vertex_sample[index]][vec[index]])
                    index = index + 1
                    

                embeddings_to_check.append(new_embedding)

        
        for embedding in embeddings_to_check:
            new_genus = check_genus(G,embedding)

            if new_genus < best_genus:
                best_embedding  = embedding
                best_genus = new_genus
                

    print(best_genus)
    return best_embedding




def check_if_each_vertex_incident_to_deg_faces(G,fixed_embedding):
    """
    Input: A sagemath graph object G, and an embedding of G
    Output: A boolean, true if the each vertex v is incident to deg(v) faces in the embedding, False otherwise
    """
    
    Faces = G.faces(embedding=fixed_embedding)

    is_good = True

    for current_face in Faces:

        if is_good == False:
            break

        vertex_list_of_face = []

        counter = 0

        while counter < len(current_face):

            vertex_list_of_face.append(current_face[counter][1])
            counter = counter + 1

        for u in vertex_list_of_face:

            if vertex_list_of_face.count(u) > 1:
                is_good = False
                break
        
    return is_good
    
    
def number_of_polyhedral_faces(G,fixed_embedding):

    """
    Input: A sagemath graph object G, and an embedding of G
    Output: An integer which counts the number of polyhdral faces in this embedding of G
    """

    Faces = G.faces(embedding=fixed_embedding)

    number_poly_faces = 0

    for current_face in Faces:

        is_good = True
        vertex_list_of_face = []

        counter = 0

        while counter < len(current_face):

            vertex_list_of_face.append(current_face[counter][1])
            counter = counter + 1

        for u in vertex_list_of_face:

            if vertex_list_of_face.count(u) > 1:
                is_good = False
                break
                
        if is_good == True:
            number_poly_faces = number_poly_faces + 1
        
    return number_poly_faces

def search_for_each_vertex_incident_to_deg_faces(G,generations, k,current_generation_embbeding):
    """
    Input: A sagemath graph object G, and an embedding of G 
    Output: Runs a genetic algorithm searching for a polyhedral embedding, optimising over k randomly sampled "columns", 
           with objective function the number of polyhedra faces, the search runs n iterations/generations. Returns the best embedding found.
    """

############# Be careful, this function runs in more than (maxDegree!)^k time
    
    best_embedding = current_generation_embbeding 
    counter = 0
    
    stopper = False

    polyhedral_found = []

    while counter < generations:

        #if stopper == True:
            #break

        counter = counter + 1
        Samples = []
        
        for i in range(10):
        
            sample  = random.sample(G.vertices(), k)
            Samples.append(sample)

        for vertex_sample in Samples:
            
            optimisation_dict = {}
            vector_list = []
            embeddings_to_check = []

            for u in vertex_sample:
                Options = list(it.permutations(G.neighbors(u)))
                optimisation_dict[u] = list(Options)
                vector_list.append(list(range(math.factorial(len(G.neighbors(u))))))

            vectors = list(it.product(*vector_list))
                

            for vec in vectors:
                
                new_embedding = copy.copy(best_embedding)
                index = 0

                while index < k:

                    # Bug here
                    new_embedding[vertex_sample[index]] = list(optimisation_dict[vertex_sample[index]][vec[index]])
                    index = index + 1
                    

                embeddings_to_check.append(new_embedding)

        
        for embedding in embeddings_to_check:

            if number_of_polyhedral_faces(G,embedding) >= number_of_polyhedral_faces(G,best_embedding):
                best_embedding  = embedding
                print([number_of_polyhedral_faces(G,best_embedding), len(list(G.faces(embedding=best_embedding)))], end='\n')

            
            if number_of_polyhedral_faces(G,embedding) == len(list(G.faces(embedding=embedding))) and embedding not in polyhedral_found:
                print("found a new polyhedral embedding \n")
                print(embedding)
                polyhedral_found.append(embedding)
                #stopper = True
                #break
                

    print([best_embedding,number_of_polyhedral_faces(G,best_embedding), len(list(G.faces(embedding=best_embedding)))])
    return polyhedral_found


def REGsearch_for_each_vertex_incident_to_deg_faces(G,generations, k,current_generation_embbeding):
    """
    Input: A sagemath graph object G, and an embedding of G 
    Output: Runs a genetic algorithm searching for a polyhedral embedding, optimising over k randomly sampled "columns", 
           with objective function the number of polyhedra faces, the search runs n iterations/generations. Returns the best embedding found.
    """

############# Slightly optimized for regualr graphs
    
    best_embedding = current_generation_embbeding 
    counter = 0
    
    D = max(G.degree())
    d = min(G.degree())

    if d!= D:
        print("Error: please give me a regular graph")
        return
        
    vectors = list(it.product(range(math.factorial(d)), repeat = k))
    
    stopper = False

    polyhedral_found = []

    while counter < generations:

        #if stopper == True:
            #break

        counter = counter + 1
        Samples = []
        
        for i in range(10):
        
            sample  = random.sample(G.vertices(), k)
            Samples.append(sample)

        for vertex_sample in Samples:
            
            optimisation_dict = {}
            embeddings_to_check = []

            for u in vertex_sample:
                Options = list(it.permutations(G.neighbors(u)))
                optimisation_dict[u] = list(Options)
                

            for vec in vectors:
                
                new_embedding = copy.copy(best_embedding)
                index = 0

                while index < k:

                    new_embedding[vertex_sample[index]] = list(optimisation_dict[vertex_sample[index]][vec[index]])
                    index = index + 1
                    

                embeddings_to_check.append(new_embedding)

        
        for embedding in embeddings_to_check:

            if number_of_polyhedral_faces(G,embedding) >= number_of_polyhedral_faces(G,best_embedding):
                best_embedding  = embedding
                print([number_of_polyhedral_faces(G,best_embedding), len(list(G.faces(embedding=best_embedding)))], end='\n')

            
            if number_of_polyhedral_faces(G,embedding) == len(list(G.faces(embedding=embedding))) and embedding not in polyhedral_found:
                print("found a new polyhedral embedding \n")
                print(embedding)
                polyhedral_found.append(embedding)
                #stopper = True
                #break
                

    print([best_embedding,number_of_polyhedral_faces(G,best_embedding), len(list(G.faces(embedding=best_embedding)))])
    return polyhedral_found