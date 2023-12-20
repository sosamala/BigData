/** @file group_1.c
 *  @brief LSH implementation
 *
 *  This contains the required implemented
 *  methods for Big Data Assigment 2 
 *
 *  @authors James Ballari, Kuljeet Kaur, Samala Samala, Tishya Sohankumar Thakkar 
 *  @bug No known bug.
 */
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<stdbool.h>
#include<float.h>
#define DEBUG 0
#define LOG 1

/** @brief Linked List Data Structure.
 *
 * used to store all the hash for data point with similar
 * hascode under a linked list, to avoid exhaustive search.
 * for a matching hash.
 *
 */
struct node {
   int data;
   struct node *next;
};


/** @brief Calculates A Random value
 *
 * Returns a random value between max and -max
 *
 * @return A random value from max to -max
 */
double getRandomValue(int max) {
    double random_value = (double) ( (double)rand() / (double)RAND_MAX ) * max;
    if(rand() & 1)  return random_value; 
    return -random_value;
}


/** @brief Calculates the euclidean distance
 *
 * Given two poionts 'pt1' and 'pt2' of dimensions 'dim'
 * this method calculates the euclidean distance.
 *
 * @return the euclidean distance
 */
double euclidean_distance(double * pt1, double * pt2, int dim){

    double sum = 0;
    for(int i = 0 ; i < dim ; i++){ 
        sum += pow( (pt1[i] - pt2[i]), 2);
    }
    return sqrt(sum);

}

/** @brief Calculates the euclidean distance
 *
 * Given two poionts 'pt1' and 'pt2' of dimensions 'dim'
 * this method calculates the euclidean distance.
 *
 * @return the euclidean distance
 */
double euclidean_distance_int(int * pt1, int * pt2, int dim){

    double sum = 0;
    for(int i = 0 ; i < dim ; i++){
        sum += pow( (pt1[i] - pt2[i]), 2);
    }
    return sqrt(sum);

}

/** @brief Calculates the hashCode for a given array.
 *
 * Given an array of size m, and a mod, calcuates a
 * hashcode 
 *
 * @return the hashcode or bucket index.
 */
int myHashCode(int m, int arr[m], int mod) {
    int result = 7;
    int prime = 31;
    for(int i=0;i<m;i++) {
        result = ( ((result * prime) %mod) + arr[i] ) % mod ;
        //printf("result = %d, arr[i] = %d\n", result, arr[i]);
    }
    if(LOG)
        //printf("m = %d, mod = %d, result = %d\n", m, mod, result);
    if(result < 0) result *= -1;
    return result;
}

/** @brief Initializes h values between -1 and 1.
 *
 * Given an array of size m x dim, populates its values
 * between -1 and 1
 *
 * @return the hashcode or bucket index.
 */
void initalize_hvalues(int dim, int m, double **h) {
    if(LOG)
        printf("\nh[%d][%d] Values:\n",m,dim);
    for(int i=0; i<m; i++) {
        for(int j=0;j<dim;j++) {
            h[i][j] = getRandomValue(1);
            if(LOG)
                printf("%0.3lf\t",h[i][j]);
        }
        if(LOG)
            printf("\n");
    }
    printf("\n");
}

/** @brief Performs LSH algorithm and creates and sorts
 * datapoints into their corresponding hash buckets.
 *
 * Given data with ndata points each of size dim and
 * Constants m, W, h, b values which are required for
 * LSH Algorithm, it calculates the hash value for each 
 * datapoint and re-arragnes the datapoints as needed.
 * popualtes the cluster_size, cluster_start and cluster_hashval
 * arrays.
 *
 * @return the total number of buckets.
 */
int LSH(int dim, int ndata, double *data, // dataset info
 int m, double W, double **h, double *b, // user-input LSH parameters
 int *cluster_start, int *cluster_size, int **cluster_hashval) {
    
    // 1. Initalize h values between -1 and 1.
    initalize_hvalues(dim, m ,h);
    
    // 2. Compute temp_hash since we do not know the number of clusters yet.
    // temp_hash[ndata] -> ith index has m points 0..m-1 which correspond to,
    // a the m hash point values. 
    // m hold total points in the cluster,
    // m+1 holds initial starting position of the cluster
    // m+2 holds the number of points assigned.
    int temp_hash_size = 0;
    int **temp_hash = (int **)malloc(ndata * sizeof(int*));
    for(int i=0;i<ndata;i++) {
        temp_hash[i] = (int *)malloc((m+3) * sizeof(int));
        temp_hash[i][m] = temp_hash[i][m+1] =
           temp_hash[i][m+2] = 0;
    }
    
    // 3. Initialize cluster_assign[ndata] stores to which cluster does
    // ith point belongs to
    int cluster_assign[ndata];
    for(int i = 0 ; i< ndata; i++){
        cluster_assign[i] = -1;
    }
    int cluster_assign_index = 0;
    
    // 4. mod - is the number of buckets created to store the hash_values
    // in a map in order to avoid exhaustive search each time.
    
    // calculating and adjusting mod values as required.
    int mod = ndata/2;
    if(mod < 2) mod = 2;
    //if( mod <2 ) mod = 2;
    //if( mod > 10000) mod = 10000;
    
    if(LOG) {
        printf("Mod value is %d\n", mod);
        printf("\n");
    }
    
    // 5. Initialize hash buckets, where the hashes with the same hashcode,
    // belong to the same hash_bucket.
    struct node *hash_buckets[mod]; 
    for(int i=0;i<mod;i++) {
       hash_buckets[i] = NULL;
    }
                            
    // 6. Perform LSH 
    for(int i = 0; i< ndata; i++) {
        
        // 6.1 Calculating the hash for the ith data point.
        int current_hash[m];
        for(int j=0; j < m; j++) {
            current_hash[j] = 0 ;
                
            for(int k = 0; k < dim ;k++) {   
                current_hash[j] += data[i*dim+k]*h[j][k];                       
            } // k
            current_hash[j] -= b[j];
            current_hash[j] /= W;
        } // j
        
        if(DEBUG) {
            printf("Current Hash for point %d is :\n",i);
            for(int l=0;l<m;l++) {
                printf("%d ",current_hash[l]);
            }
            printf("\n");
        }
        
        
        
        // 6.2 calculating the bucket index for the hash.    
        int hash_bucket_index = myHashCode(m, current_hash,mod);
        struct node *current_bucket = hash_buckets[hash_bucket_index];
        
        if(DEBUG) {
            printf("current bucket index is %d \n", hash_bucket_index);
        }
        
        if(current_bucket == NULL) { // 6.2.1 no match - first entry..
            cluster_assign[i] = temp_hash_size;
            temp_hash[temp_hash_size][m] += 1;
            hash_buckets[hash_bucket_index] = (struct node*)malloc(sizeof(struct node));
            hash_buckets[hash_bucket_index]->data = temp_hash_size;
            hash_buckets[hash_bucket_index]->next = NULL;
            for(int k = 0 ; k < m ; k++){
                    temp_hash[temp_hash_size][k] = current_hash[k];  
            }
            temp_hash_size++;
        } else { // 6.2.2 check for a matching entry..
            bool isMatch = false;
            while(true) {
                int toCheckHashIndex = current_bucket->data;
                bool hashMatched = true;
                for(int k = 0 ; k < m; k++){
                    if(temp_hash[toCheckHashIndex][k] != current_hash[k]) {
                        hashMatched = false;   
                        break;
                    }
                }
                if(hashMatched) { // 6.2.2.1 found a matching hash.
                    isMatch = true;
                    cluster_assign[i] = toCheckHashIndex;
                    temp_hash[toCheckHashIndex][m] += 1;
                    //temp_hash_size++;
                    break;
                }
                    
                if(current_bucket->next != NULL) {
                    current_bucket = current_bucket->next;
                } else {
                    break;
                }               
            }
            
            if(!isMatch) { // 6.2.3 no match - last entry - create a new node
                struct node *new_node = (struct node*)malloc(sizeof(struct node));;
                new_node->next = NULL;
                new_node->data = temp_hash_size;
                temp_hash[temp_hash_size][m] += 1;
                current_bucket->next = new_node;
                for(int k = 0 ; k < m ; k++){
                    temp_hash[temp_hash_size][k] = current_hash[k];  
                }
                cluster_assign[i] = temp_hash_size;
                temp_hash_size++;
                
            }
        
        } 
    } // 
    
    
    if(LOG) {
        printf("Printing Hash Buckets\n");
        for(int i=0;i<mod;i++) {
            printf("%d. ", (i+1));
            struct node *current_bucket = hash_buckets[i];
            if(current_bucket != NULL) {
                printf("%d\t", current_bucket->data);
                while(current_bucket->next != NULL ) {
                    current_bucket = current_bucket->next;
                    printf("%d\t", current_bucket->data);
                }
            } else {
                printf("Bucket is empty");
            }
            printf("\n");
        }
        printf("\n");
    }
    
    // 7. Calculating temp_hash[..][m+1] values..
    for(int i=1;i<temp_hash_size;i++) {
        temp_hash[i][m+1] = temp_hash[i-1][m] + temp_hash[i-1][m+1];
    }
    
   
    
    // 8. move points according to their Bucket position 
    double tempPoints[ndata*dim];
    for(int i=0;i<ndata;i++) {
        int cluster = cluster_assign[i];
        int position = temp_hash[cluster][m+1]+temp_hash[cluster][m+2];
        for(int j=0;j<dim;j++) {
            tempPoints[position*dim+j] = data[dim*i+j];
        }
        temp_hash[cluster][m+2] += 1;
        
    }
    for(int i=0;i<ndata*dim;i++) {
            data[i] = tempPoints[i];
    }
    
     

    if(LOG){
        printf("Print temp hash Assign\n");
         for(int i = 0; i < temp_hash_size; i++){
             printf("%d.\t",i);
             for(int j = 0 ; j < m+3; j++) {
              printf(" %d\t",   temp_hash[i][j]);
             }
             printf("\n");
         
        }
        
        printf("\nPrint cluster Assign\n");
        for(int i = 0;i < ndata;i++){
            printf(" %d",cluster_assign[i]);   
        }
        printf("\n\n");
    }
    
    if(DEBUG) {
        printf("Printing All the points\n");
        for(int i=0;i<ndata;i++) {
            printf("%d. ", i+1);
            for(int j=0;j<dim;j++) {
                printf("%0.2lf, ", data[dim*i+j]);
            }
            printf("\n");
        }
        printf("\n");
    }
    
    // 9. Set cluster size and cluster size arrays.
    for(int i=0;i<temp_hash_size;i++) {
        cluster_start[i] = temp_hash[i][m+1];
        cluster_size[i] = temp_hash[i][m+2];
    }
    
    if(LOG) {
        printf("Cluster start and Cluster Size is\n");
        for(int i=0;i<temp_hash_size; i++) {
                printf("(%d, %d), ", cluster_start[i], cluster_size[i]);
        }
        printf("\n");
    }
    
    // 10. Finally set the cluster_hashval values.
    for(int i=0;i<temp_hash_size;i++){
        cluster_hashval[i] = (int *) malloc(m * sizeof(int));
        for(int j=0;j<m;j++) {
                cluster_hashval[i][j] = temp_hash[i][j];
            }
    }
    
    // return the total number of clusters.
    return temp_hash_size;
}


/** @brief Perform exhaustive search on a cluster.
 *
 * Given a 'cluster' index,  'cluster_start[cluster]' 
 * and 'cluster_size[cluster]' indexing the 
 * points in 'data' points with 'dim' dimensions, performs an exhaustive 
 * search of all points to find the minimum distance points.
 *
 * @return the minimum distance and update the result_pt
 */
double search_cluster(
    int cluster, int dim, double *data, int * cluster_start,
    int * cluster_size,double 
    *query_pt, double *result_pt ){
    
        double nearest_distance = DBL_MAX;
        int current_cluster_size = cluster_size[cluster];
        int current_cluster_start = cluster_start[cluster];
    
        for(int i = 0; i < current_cluster_size; i++) {

            double * curr_pt 
                =  &data[(current_cluster_start+i)*dim];
            double distance 
                = euclidean_distance(query_pt, curr_pt, dim);
                
            if(DEBUG) {
                printf("\nDistance: %lf  %lf\n", distance,
                nearest_distance);
            }
            //nearest_distance = min(nearest_distance, distance);
            if(distance < nearest_distance){
                nearest_distance = distance;
                for(int l=0; l < dim ;l++){
                    result_pt[l] = curr_pt[l];
                    if(DEBUG) {
                        printf("%lf ", curr_pt[l]);
                    }
                }
                if(DEBUG) {printf("\n");}
                
            }
        }
        if(LOG) {
                printf("\nNearest Distance: %lf\n",nearest_distance);
        }
        return nearest_distance;
}


/** @brief Searches the clusters using LSH algorithm.
 *
 * Calcualte the hash value for query point.
 * search for the closest cluster, the cluster with a 
 * hash value that corresponds to the minimum 
 * distance form the query hash. 
 * exhasutively search all points in the cluster for the
 * nearest point. sets the result_pt
 *
 * @return 0.
 */
int search_LSH(int dim, int ndata, double *data,
 int m, double W, double **h, double *b, 
 int nclusters, int *cluster_start, int *cluster_size, 
 int **cluster_hashval, double *query_pt, double *result_pt) {

    // 1.Calculate hash for query point.
    int current_hash[m];
    for(int j=0; j < m; j++) {
        current_hash[j] = 0 ;   
        for(int k = 0; k < dim ;k++) {   
            current_hash[j] += query_pt[k]*h[j][k];                       
        } // k
        current_hash[j] -= b[j];
        current_hash[j] /= W;
        //printf("\n Current Hash %d : %d", j,current_hash[j]);
    } // j
    
    if(LOG) {
        printf("\nQuery point hashvalue is :");
        for(int i=0;i<m;i++) {
            printf("%d ",current_hash[i]);
        }
        printf("\n\n");
    }
    
    // 2. go through all the clusters, get the min euclidean distance cluster.
    // if min distance = 0 break;
    int choosen_cluster[ndata];
    int choosen_cluster_size = 0;
    double min_distance = DBL_MAX;
    for(int i=0;i<nclusters;i++) {
        double distance = euclidean_distance_int(
            cluster_hashval[i],current_hash,m);
        if(LOG) {
            printf("Distance from the %d Cluster is : %0.4lf\n",i, distance );
        }
        if(DEBUG){
            for(int t = 0; t < dim;t++){
            printf("\nCluster %d Cluster_hashval %d, current_hash %d",i,cluster_hashval[i][t],current_hash[t]);
            }
        }
        if(distance < min_distance) {
            
            if(DEBUG){
             printf("\n Distance: Cluster %d distance %d",i,distance);   
            }
            min_distance = distance;
            choosen_cluster_size = 1;
            choosen_cluster[0] = i;
            if(LOG)
                printf("New Minimum Distance is : %0.4lf, selected Cluster is %d\n", 
                       min_distance, i);
            if(min_distance == 0) {
                break;
            }
        } else if(distance == min_distance){
            choosen_cluster[choosen_cluster_size++] = i; 
        }
        
        
    }
    
    double nearest_cluster_dist = DBL_MAX;
    double ptr[dim] ;
    // exhaustively search all the points in the cluster.
    for(int i = 0 ; i < choosen_cluster_size; i++){
        
    double cluster_dist = search_cluster(choosen_cluster[i], dim, data, cluster_start,
     cluster_size, query_pt, ptr);
        printf("\nExhaustively Searching cluster %d which has the closest distance %lf\n",choosen_cluster[i], cluster_dist);
    
        if(cluster_dist <  nearest_cluster_dist){
            for(int i = 0; i < dim ; i++){
                    result_pt[i] = ptr[i];
                }
            nearest_cluster_dist = cluster_dist;
        }
    }
    
    return 0;// return the result point
    
    
}



void test_getRandomValue(){
    for(int i=0;i<20;i++) 
        printf("%lf, ", 
               (double)getRandomValue(1));
}

/** @brief The main method.
 *
 * calls the LSH method and search_LSh method with
 * the required data.
 *
 * @return 0
 */
int main() {
    //srand(time(NULL); // randomize seed.
    srand(0); // same seed.
    int dim = 4;
    int ndata = 32;
    double data[] = { 
        1,2,13,4,
        0,-2,4,5,
        5,2,1,3,
        -3,-1,4,1,
        -13,2,2,3,
        
        -4,3,-1,2,
        12,2,3,4,
        3,4,1,-1,
        2,3,1,-1,
        -1,2,3,4,
        
        2,-3,2,1,
        4,1,-3,2,
        1,-3,1,-4,
        -4,3,2,5,
        2,3,4,3,
        
        -2,5,3,3,
        -4,2,-5,4,
        5,5,1,3,
        1,1,0,1,
        4,1,5,1,
        
        5,3,-4,-1,
        -4,-1,-2,-3,
        4,2,-3,3,
        0,-5,1,-1,
        4,-4,-2,4,
        
        0,1,15,5,
        -4,0,-1,-2,
        4,-3,-3,1,
        -1,2,4,-3,
        12,1,-5,2,
        
        5,-33,-1,-3,
        4,-3,-1,5
    };
    
    double query_pt[] = {5,1,-15,-1};
    
    int m = 3; // set value for m.
    double W = 3;
    double **h;
    // Initilaize h,
    h = (double **) malloc(m * sizeof(double *));
    for(int i=0; i<m; i++) {
        h[i] = (double *)malloc(dim * sizeof(double));
    }
    
    double b[] = {1.4,0.2,-0.4};
    int cluster_start[ndata];
    int cluster_size[ndata];
    int **cluster_hashval;
    cluster_hashval = (int **) malloc(ndata * sizeof(int *));
    
    
    
    
    int nclusters = LSH(dim, ndata, data, m, W, h, b, 
        cluster_start, cluster_size, 
        cluster_hashval);
    
    printf("\nNumber of Clusters = %d \n", nclusters);
    

    
    //double query_pt[] = {3,2,4};
    //double query_pt[] = {1,2,3};
    //double query_pt[] = {0,-2,5};
    
    double result_pt[dim]; 
    if(LOG) {
        printf("\nQuery Point\n");
        for(int i = 0; i<  dim ; i++){
            printf("%0.3lf ",query_pt[i]);
        }
        printf("\n");
    }
    
    search_LSH( dim,  ndata, data, m,  W, h, b, 
      nclusters, cluster_start, cluster_size, 
     cluster_hashval, query_pt, result_pt);
    
    printf("\n\nResult! \n");
    for(int i = 0; i<  dim ; i++){
        printf("%0.3lf ",result_pt[i]);
    }
    
    
    printf("\n\nFin! \n");
    
    return 0;
}
