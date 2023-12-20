/** @file group_1.c
 *  @brief Kmeans implementation
 *
 *  This contains the required implemented
 *  methods for Big Data Assigment 3 
 *
 *  @authors James Ballari, Kuljeet Kaur, Samala Samala, Tishya Sohankumar Thakkar 
 *  @bug No known bug.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<float.h>
#include<time.h>
#include<stdbool.h>
#define DEBUG 0
#define LOG 0



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

/** @brief returns k initial centers for the K-means
 *
 * Selects k centers which are equally far away from 
 * each other.
 *
 * @return the number of centers
 */
int initial_centers(int dim, int i0, int im, double *data, int k, 
                    double *cluster_centroid) {
    
    // 1. Select Cluster Center points from datum
    // 1.1 Choose inital cluster center.
    int ndata = im - i0 + 1;
    int cluster_center_points[k];
    bool *cluster_selected_points = (bool*) malloc(ndata * sizeof(bool));
    for(int i=0;i<ndata;i++) cluster_selected_points[i] = false;
    
    cluster_center_points[0] = (int)(rand() % ndata);
    cluster_selected_points[cluster_center_points[0]] = true;
    
    if(DEBUG) printf("First point is %d \n",cluster_center_points[0] );
    if(k == 1) return k;
    
    // 1.2 choosing farthest point from 1st selected cluster.
    double max_distance = -1;
    for(int i=0;i<ndata;i++) {
        if(!cluster_selected_points[i]) {
            double distance 
                = euclidean_distance(
                    &data[cluster_center_points[0]*dim],
                    &data[i*dim],
                    dim);
            if(distance > max_distance) {
                max_distance = distance;
                cluster_center_points[1] = i;
                if(DEBUG) printf("max distance  = %lf\n", distance);
            }
        }
    }
    cluster_selected_points[cluster_center_points[1]] = true;
    if(DEBUG) printf("Second point is %d \n",cluster_center_points[1] );
    if(k == 2) return k;
    
    // 1.3 choosing the rest of the cluster centers 
    printf("..");
    int j =2;
    while(j<k) {
        max_distance = -1;
        for(int i=0;i<ndata;i++) {
            if(!cluster_selected_points[i]) {
                double min_distance = DBL_MAX;
                for(int k=0; k<j; k++) {
                    double distance = euclidean_distance(
                    &data[cluster_center_points[k]*dim],
                    &data[i*dim],
                    dim);
                    if(distance < min_distance) {
                        min_distance = distance;
                    }
                }
                if(min_distance > max_distance) {
                    cluster_center_points[j] = i;
                    max_distance = min_distance;
                }
            }
        }
        cluster_selected_points[cluster_center_points[j]] = true;
        j++;
        printf(".");
    }
    printf("\n");
    
    // 1.4 updating the cluster centeriod accordingly.
    for(int i=0;i<k;i++) {
        int centroid_point_index = cluster_center_points[i];
        for(int j=0;j<dim;j++) {
            cluster_centroid[i*dim+j] = data[centroid_point_index * dim + j];
        }
    }
    
    // 1.5 Printing all the selected cluster centers.
    if(LOG) {
        printf("\nSlected Cluster center points are : \n");
        for(int i=0;i<k;i++) {
            printf("%d. %d -> ", i+1, cluster_center_points[i]);
            for(int j=0;j<dim;j++) {
                
                
                printf("%0.0lf, ", cluster_centroid[i*dim+j]);
            }
            if(LOG) printf("\n");
        }
    }
    
    return k;
    
}



/** @brief Cretes the KMeans Clusters
 *
 * Creates KMeans Clusters for 'ndata' number of points
 * with 'dim' dimensions whose data is found in 'data'.
 * creates 'k' no. of clusters, where each cluster is of 
 * 'cluster_size[i]' size starting from 'cluster_start[i]'
 * position in 'cluster_assign', where the centroid of each
 * cluster is in 'cluster_centroid' and 
 * cluster radius is in 'cluster_radius'.
 *
 * @return no return values.
 */
/********************************************************************
 Array Sizes:
 cluster_centroid[k*dim]: input -- stores initial k cluster centers
                          output-- stores final k cluster centers

 cluster_radius[k]:output
                   Stores the radius of each output cluster

 cluster_start[k]: output
                   Stores the index of each cluster's starting datum

 cluster_size[k]:  output
                   Stores the num of data points in each cluster

 cluster_assign[ndata]: buffer that stores the membership of data 
*********************************************************************/
int kmeans(int dim, int i0, int im, double *data, int k, 
           int *cluster_start, int *cluster_size,
           double *cluster_radius, double *cluster_centroid,
           short *cluster_assign)   {
    
    int ndata = im-i0+1;
    // 1. initial_centers
    k = initial_centers(dim, i0, im, data, k, 
                    cluster_centroid);
    
    // 2. K means Algorithm.
    bool stop_iteration = false;
    int max_tries = ndata*k, tries = 0;
    if(LOG) printf("\nK means clustering Starts\n");
    while(!stop_iteration && tries  < max_tries) {
        
        int count_cluster_change = 0;
        if(LOG) printf("*******************************************");
        if(LOG) printf("\nTries = %d\n", tries+1);
        // setting cluster size to zero.
        for( int i=0;i<k;i++) {
            cluster_size[i] = 0;
        }
        if(LOG) printf("\nChanges = ");
        for(int i=0;i<ndata;i++) {
            int old_cluster = cluster_assign[i];
            int new_cluster = -2;
            double min_distance = DBL_MAX;
            
            for(int j=0;j<k;j++) {
                double distance = euclidean_distance(
                    &cluster_centroid[j*dim],
                    &data[i*dim],
                    dim);
                if(distance < min_distance) {
                    min_distance = distance;
                    new_cluster = j;
                }
            }
            if(old_cluster != new_cluster) {
                count_cluster_change++;
                if(LOG) printf("[%d:%d->%d],",i,old_cluster,new_cluster);
                
            }
            cluster_assign[i] = new_cluster;
            cluster_size[new_cluster]++;
        }
        
        // calculate centroid of each cluster.
        for(int i=0;i<k;i++) {
            for(int j=0;j<dim;j++) {
                cluster_centroid[i*dim+j] = 0;
            }
        }
        if(LOG) printf("\n\nCluster Assign = ");
        for(int i=0;i<ndata;i++) {
            if(LOG) {
            printf("(%d -> %d),",i, cluster_assign[i]);
            }
            for(int j=0;j<dim;j++) {
                cluster_centroid[cluster_assign[i]*dim+j] += 
                    data[i*dim + j];
            }
        }
        if(LOG) printf("\n");
        
        
        for(int i=0;i<k;i++) {
            for(int j=0;j<dim;j++) {
                cluster_centroid[i*dim+j] /= cluster_size[i];
            }
        }
        
        if(LOG) {
            printf("\nCluster Centroids: \n");
            for(int i=0;i<k;i++) {
                printf("%d. ", i+1);
                for(int j=0;j<dim;j++) {
                    printf("%0.4lf, " , cluster_centroid[i*dim+j] );
                }
                printf("\n");
            }
        
        }
        
        // break conditions
        if(count_cluster_change == 0) stop_iteration = true;
        tries++;
        printf(".");
    }
    printf("\n");
    // 3. Populate cluster start, re-arrage data, find radius.
    
    // 3.1 Populate cluster start from cluster size.
    cluster_start[0] = 0;
    for(int i=1;i<k;i++) {
        cluster_start[i] = cluster_start[i-1]+cluster_size[i-1];
    }
    
    if(LOG) {
        printf("\nCluster start and Cluster Size is\n");
        for(int i=0;i<k;i++) {
             printf("(%d, %d), ", cluster_start[i], cluster_size[i]);
        }
        printf("\n");
    }
    
    
    
    // 3.2 move points according to their cluster position
    double *tempPoints = (double*) malloc(ndata*dim * sizeof(double));
    double cluster_fill_size[k];
    memset(cluster_fill_size,0, sizeof(cluster_fill_size));
    for(int i=0;i<ndata;i++) {
        int cluster = cluster_assign[i];
        int position = cluster_start[cluster]+cluster_fill_size[cluster];
        for(int j=0;j<dim;j++) {
            tempPoints[position*dim+j] = data[dim*i+j];
        }
        cluster_fill_size[cluster]+= 1;
        
    }

    if(DEBUG) printf("\n Points are: \n");
    for(int i=0;i<ndata*dim;i++) {
            data[i] = tempPoints[i];
            if(DEBUG) {
                if(i%dim ==0) if(LOG) printf("\n");
                printf("%0.2lf, ",data[i]);
                }
    }
    
    // 3.3 Calculate Cluster Radius.
    // Distance from the Centroid to the Farthest Point in the Cluster
    for (int i=0;i<k;i++) {
        int clu_start = cluster_start[i];
        int clu_end = cluster_start[i] + cluster_size[i] - 1;
        cluster_radius[i] = 0;
        for(int j=clu_start;j<=clu_end;j++) {
            double distance = euclidean_distance(
                &cluster_centroid[i*dim], 
                &data[j*dim], 
                dim);
            if(distance > cluster_radius[i]) {
                cluster_radius[i] = distance;
            }
        } 
    }
    
    if(LOG) {
        printf("\nCluster Radius:\n");
        for(int i=0;i<k;i++) 
            printf("%d) %lf\n",i+1, cluster_radius[i]);
    }
    
    return k;
    
}


/** @brief Calculates the minimum
 *
 * Given two values 'a' and 'b'
 *
 * @return the minimum of a and b
 */
double min(double a, double b){
    return a<b ? a : b;
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
    
        for(int i = 0; i < current_cluster_size; i++){

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
        if(DEBUG) {
                printf("\nNearest Distance: %lf\n", nearest_distance);
        }
        return nearest_distance;
}

/** @brief search for the nearest point using clusters of KMeans.
 *
 * Given 'ndata' number of data points, each with 'dim' dimensions, 
 * search for 'query_pt', using 'k' clusters of KMeans with
 * corresponding 'cluster_start', 'cluster_size', 'cluster_centroid'
 * and cluster_radius.
 * 
 * The Search has two steps. 
 *   1. find the cluster with the closest centroid, perform exhaustive
 *      search to find the closet point in the centroid, call it dmin.
 *   2. if dmin ==0 we have the search point else
 *        for each cluster 
 *          calcuate if dmin > (distance to centroid - Radius of the cluster)
 *          then perform exhaustive serach in that cluster calcuate new_d_min.
 *           if new_d_min is less than dmin, update dmin = new_d_min. 
 *                                          and update result point.
 * @return no return value, updates result_pt. */
/***************************************************************
 * search_kmeans() returns the number of data points checked   *
 * query_pt[dim]:  input, the query point                      *
 * result_pt[dim]: output, the closest data point to the query *
 ***************************************************************/

int search_kmeans(int dim, int ndata, double *data, int k,
           int *cluster_start, int *cluster_size, 
           double *cluster_radius, double *cluster_centroid,
           double *query_pt, double *result_pt) {
    
    int is_result = 0;
    double nearest_distance = DBL_MAX ;
    double nearest_pt[dim];

        
    // search in the cluster with the closest centroid to query point.

    double nearest_ctr_dist = DBL_MAX ;
    int nearest_ctr[k];
    int nearest_ctr_size = 0;
    bool checkedclusters[k];
    for(int i=0;i<k;i++) checkedclusters[i] = false;
    
    //int nearest_ctr = -1;
    double dist_qp_ctr[k];
    for(int cltr = 0; cltr < k; cltr++) {
        double *centroid = &cluster_centroid[cltr*dim];
        dist_qp_ctr[cltr] = euclidean_distance(query_pt, centroid, dim);
        if(DEBUG) {
            printf("\nQuery Point\n");
            for(int i = 0; i<  dim ; i++){
                printf("%lf, ",centroid[i]);
            }
        }
        //nearest_ctr_dist = min(distance,nearest_ctr_dist);
        if(dist_qp_ctr[cltr] < nearest_ctr_dist ) {
            nearest_ctr_size = 1;
            nearest_ctr[0] = cltr;
            nearest_ctr_dist = dist_qp_ctr[cltr];
        } else if(dist_qp_ctr[cltr] == nearest_ctr_dist ) {
             nearest_ctr[nearest_ctr_size++] = cltr;
        } 
    }

    if(LOG) {
        printf("\nNearest Cluster Distance is : %lf,\nNearest Clusters are: ", nearest_ctr_dist);
        for(int i=0;i<nearest_ctr_size;i++){
            printf("%d, ", nearest_ctr[i]);
        }
        printf("\n");
    }
    
    double nearest_distance_temp = DBL_MAX;
    double d_min;
    int nearest_cluster_index = 0;
    for(int i=0;i<nearest_ctr_size;i++){
        if(LOG) printf("\nSearching in Cluster = %d for the closest point.\n", nearest_ctr[i]);
         checkedclusters[nearest_ctr[i]] = true;
         d_min = search_cluster(nearest_ctr[i], dim, data,cluster_start,
            cluster_size,query_pt,nearest_pt);
        if(LOG) printf("Nearest Point is : ");
        for(int j=0;j<dim;j++){
            if(LOG) printf("%lf, ", nearest_pt[j] );
        }
        if(LOG) printf("\nDistance form Query point is %lf\n", d_min);
        if(d_min < nearest_distance_temp) {
            nearest_distance_temp = d_min;
            for(int j = 0; j < dim ; j++){
                    result_pt[j] = nearest_pt[j];
            }
            nearest_cluster_index = i;
            
        }
    }
    
    d_min = nearest_distance_temp;
    
    if(LOG) printf("\nD_min %lf\n",d_min); 
    if(d_min == 0.0) 
        return 0;
    if(LOG) printf("\nChecking Other clusters for possible nearest point.\n");
    
    for(int cltr = 0; cltr < k; cltr++) {
        if(LOG) printf("*************************************************************");
        if(checkedclusters[cltr]) {
            if(LOG) printf("\nSkipping Cluster %d, as its been already searched.\n", cltr);
            continue;
        }
        if(LOG) printf("\nCluster=%d \n\t Cluster Radius=%lf, "
               "\n\t Distance from Query point to Cluster Centroid = %lf"
               "\n\t Current Minimum Distance = %lf",
               cltr,cluster_radius[cltr],dist_qp_ctr[cltr], d_min ); 
        if(d_min > (dist_qp_ctr[cltr] - cluster_radius[cltr])){
            if(LOG) printf("\n\n\t Searching cluster %d for possible nearest point",cltr);
            double new_d_min = search_cluster(cltr, dim, data,cluster_start,
            cluster_size,query_pt,nearest_pt);
            if(LOG) printf("\n\t Closest Point is %lf away", new_d_min);
            if(LOG) printf("\n\t Closest Point is : ");
            for(int j=0;j<dim;j++){
                if(LOG) printf("%lf, ", nearest_pt[j] );
            }
            if(LOG)  printf("\n");
            if(new_d_min < d_min){
                for(int i = 0; i < dim ; i++){
                    result_pt[i] = nearest_pt[i];
                }
                d_min = new_d_min;
                if(LOG)  printf("\nNew D_min value %lf in cluster %d",d_min,cltr);
            }          
        } else {
            if(LOG)  printf("\n\t Unlikely for a Nearest Point to be in this cluster.");
        }
        if(LOG)  printf("\n");
          
    }
    
    
    return is_result;
}


/** @brief The main method.
 *
 * calls the KMeans method and search_kmeans method with
 * the required data.
 *
 * @return 0
 */
int main() {

    srand(0); // same seed.
    int dim = 8;
    int ndata = 1000000;
    //int ndata = 1000;
    int nqueries = 1000;
    int i0 = 0 ;
    int im = ndata-1;
    double *data = (double*) malloc(ndata*dim * sizeof(double));
    double *query_points = (double*) malloc(nqueries*dim * sizeof(double)); 
    
    int k = 1000;
    printf("Starting Generating %d number of data points  k = %d\n", ndata,k);
    
    for(int i=0;i<ndata*dim;i++) {
        data[i] = rand()%10;
    }
    

    for(int i=0;i<nqueries*dim;i++) {
        query_points[i] = rand()%10;
    }
    
    printf("Generated %d number of data points and %d number of query points, k = %d\n", ndata, nqueries,k);
    

    int cluster_start[k];
    int cluster_size[k];
    double cluster_radius[k];
    double cluster_centroid[k*dim];
    short *cluster_assign = (short*) malloc(ndata * sizeof(double));
    
    // defaulting all the values
    memset(cluster_start,-1, sizeof(cluster_start));
    memset(cluster_size,-1, sizeof(cluster_size));
    memset(cluster_radius,0, sizeof(cluster_radius));
    memset(cluster_centroid,0, sizeof(cluster_centroid));
    memset(cluster_assign,-1, sizeof(cluster_assign));
    
    
    double time_spent = 0.0;
    clock_t begin = clock();
    
    // Call Kmeans method.
    printf("Started k means\n");
    clock_t begin_kmeans = clock();
    
    kmeans(dim, i0, im, data, k,
           cluster_start, cluster_size,
           cluster_radius, cluster_centroid,
           cluster_assign);
    clock_t end_kmeans = clock();
    double keamns_time_spent = (double)(end_kmeans - begin_kmeans) / CLOCKS_PER_SEC;
    printf("\nThe elapsed time to search is %f seconds\n", keamns_time_spent);
    
    printf("End k means\n");


    //double query_pt[] = {3,2,4};
    //double query_pt[] = {-4,3,-1,2};
    //double query_pt[] = {50,-1,5,-23};
    double query_pt[dim];
    //double query_pt[] = {5,2,1};
    double result_pt[dim];
    
    if(LOG) {
        printf("\nQuery Point\n");
        for(int i = 0; i<  dim ; i++){
            printf("%lf ",query_pt[i]);
        }
    }
    if(LOG) printf("\n");
    
        
    //call search kdtree method.
    clock_t begin_search = clock();
    for(int i=0;i<nqueries;i++) {
    
        for(int j=0;j<dim;j++)
            query_pt[j] = query_points[i*dim+j];
        
        int a = search_kmeans(dim, ndata, data, k,
          cluster_start, cluster_size, cluster_radius,cluster_centroid,
          query_pt, result_pt);
        //if(i%100==0) {
        //    printf("Processed %d number of queries.\n",i);
        //}
        for(int k = 0; k<  dim ; k++)
            printf("%lf ",result_pt[k]);
        printf("\n");
    }

    clock_t end_search = clock();
    double search_time_spent = (double)(end_search - begin_search) / CLOCKS_PER_SEC;
    printf("\nThe elapsed time to search is %f seconds\n", search_time_spent);
    
   
    if(LOG)  {
    printf("\n\nResult! \n");
    for(int i = 0; i<  dim ; i++){
        printf("%lf ",result_pt[i]);
    }
    }
    
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nThe elapsed time is %f seconds\n", time_spent);
    printf("\n\nFin! \n");

    return 0;
}

// To run this program
//gcc group_1.c -o kd -lm -w && ./kd
