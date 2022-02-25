#include "Utilities.h"

// install CUnit in Ubuntu: sudo apt-get install libcunit1 libcunit1-doc libcunit1-dev
#include <CUnit/CUnit.h>   
#include "CUnit/Basic.h"

using namespace std;


int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }


// function definition for test continuous frechet between 2 curves
void test_Continuous_Frechet_Distance();  
void test_Discrete_Frechet_Distance();


int main()
{
    CU_pSuite pSuite1, pSuite2 = NULL;

    // Initialize CUnit test registry
    if (CUE_SUCCESS != CU_initialize_registry())
    return CU_get_error();

    // Add suite1 to registry
    pSuite1 = CU_add_suite("Test Continuous Frechet Distance", init_suite, clean_suite);
    if (NULL == pSuite1) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    // add test_Continuous_Frechet_Distance to suite1
    if ((NULL == CU_add_test(pSuite1, "\n\n……… Testing Continuous_Frechet_distance function……..\n\n", test_Continuous_Frechet_Distance))) {
        CU_cleanup_registry();
        return CU_get_error();
    }


    pSuite1 = CU_add_suite("Test Discrete Frechet Distance", init_suite, clean_suite);
    if (NULL == pSuite1) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    // add test Discrete_Frechet_distance to suite1
    if ((NULL == CU_add_test(pSuite1, "\n\n……… Testing Discrete_Frechet_distance function……..\n\n", test_Discrete_Frechet_Distance))) {
        CU_cleanup_registry();
        return CU_get_error();
    }



    CU_basic_run_tests();   // OUTPUT to the screen
    CU_cleanup_registry();  //Cleaning the Registry
    return CU_get_error();

}





// function declaration for test
void test_Continuous_Frechet_Distance() {
    double distance = 0.0; 

    vector<double> curve1 {10, 10};
    vector<double> curve2 {20, 20};
    distance = Continuous_Frechet_distance(&curve1, &curve2);   // output should be 10, between 2 parallel lines with distance 10
    CU_ASSERT(10 == distance); 


    vector<double> curve3 {1.5, 1.5, 1.5};
    vector<double> curve4 {2.5, 2.5};
    distance = Continuous_Frechet_distance(&curve3, &curve4);
    CU_ASSERT(1 == distance);  
}




// function declaration for test
void test_Discrete_Frechet_Distance() {
    double distance = 0.0; 

    vector<double> curve1 {10, 10};
    vector<double> curve2 {20, 20};
    vector<double> mean_curve;
    distance = Discrete_Frechet_distance(&curve1, &curve2, mean_curve);   // output should be 10, between 2 parallel lines with distance 10
    CU_ASSERT(10 == distance); 


    vector<double> curve3 {1.5, 1.5, 1.5};
    vector<double> curve4 {2.5, 2.5};
    distance = Discrete_Frechet_distance(&curve3, &curve4, mean_curve);
    CU_ASSERT(1 == distance);  
}


