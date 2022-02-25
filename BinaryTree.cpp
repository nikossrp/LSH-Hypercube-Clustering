#include "BinaryTree.h"
#include <queue>

using namespace std;

template <typename T>
BinaryTree<T>::BinaryTree(int num_curves)
{
    this->root = NULL; 
    this->height = ceil(log2(num_curves));
    this->num_nodes = pow(2, this->height)*2 - 1;       // we must create 2 ^ (num_nodes -1) - 1 nodes, we need to have the curves at the leafs only  
    int n_nodes = this->num_nodes;
    // cout << "Num of nodes: " << n_nodes << endl;

    // create an empty binary treem we will use this to find the mean curve in a cluster
    while(n_nodes--) {
        TreeNode* new_node = (TreeNode*) malloc(sizeof(TreeNode));
        new_node->right = NULL;
        new_node->left = NULL;
        new_node->parent = NULL;
        new_node->index_node = n_nodes;
        new_node->curve = NULL;
        this->insert(new_node);
    }

    // initialize leafs, we will insert curves only in leafs, other nodes is for mean curves
    init_leafs(this->root);
}




// Delete the binary tree
template <typename T>
void BinaryTree<T>::destroy_tree(TreeNode* node)
{

    if(node !=NULL)
    {
        destroy_tree(node->left);
        destroy_tree(node->right);

        // cout << "Deleting node: " << node->index_node << endl;
        if (node->curve != NULL) {
            // cout << "Curve: ";
            // for (auto c : *node->curve)
            //     cout << c << " ";
            // cout << endl;
            delete node->curve;
        }
        free(node);

    }
}
template <typename T>
BinaryTree<T>::~BinaryTree()
{
    TreeNode* root = this->get_root();
    this->destroy_tree(root);
}





template <typename T>
void BinaryTree<T>::get_mean_curve(vector<T>& mean_curve)
{
    vector<T> *curve1, *curve2;
    TreeNode* curr_node = NULL;
    int num_iterations = this->get_number_of_leafs();
    int num_iterations_temp = num_iterations;
    // compute the mean curve between 2 leafs and insert in the parent node

    vector<TreeNode*> inner_leafs;
    // this->print2DUtil(this->get_root(), 0);

    // cout << "I am inside to mean curve function\n";
    
    while(num_iterations > 1) {  // if leafs == 1 we are in the root and we have done
        
        
        // cout << "Number of leafs: " << num_iterations_temp << endl;

        for (int i = 0; i < num_iterations_temp; i+=2) {

            if (num_iterations_temp % 2 != 0) {  // if number of leafs is odd we will insert the last element to the parent instantly
                num_iterations_temp--;
            }

            curr_node = this->leafs.at(i);
            // cout << "Current leaf: " << curr_node->index_node << " When number of leafs: " << num_iterations << endl;
            // this->print_leafs();

            // cout << " (" << num_iterations_temp << "/" << num_iterations << ") We calculate mean curve between : i:" << i << " : " << this->leafs.at(i)->index_node <<  " with i:" << i+1 << " " << this->leafs.at(i+1)->index_node << endl;
            curve1 = this->leafs.at(i)->curve;



            // cout << "i: " << i << " i+1: " << i+1 << endl;
            // else we calculate the mean between 2 curves 
            curve2 = this->leafs.at(i+1)->curve;

            vector<T> mean_curve_temp;
            mean_curve_temp.push_back(0);
            Discrete_Frechet_distance(curve1, curve2, mean_curve_temp);
            // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
            // cout << "Mean curve: size: " << mean_curve_temp.size() << " cods:" << mean_curve_temp.at(0) << " " << mean_curve_temp.at(1) << endl;
            
            curr_node->parent->curve = new vector<T>();

            if (curr_node->parent->index_node != this->leafs.at(i+1)->parent->index_node) {
                cout << "Something gone wrong with get_mean_curve method\n";    // debugging
                exit(EXIT_FAILURE);
            }
            // cout << "~~~ Inserted to " << curr_node->parent->index_node << "\n";
                
            for (auto c : mean_curve_temp)
                curr_node->parent->curve->push_back(c);

            // for (auto c : *curr_node->parent->curve)
            //     cout << c << " ";
            // cout << " == ";
            // for (auto c : mean_curve_temp)
            //     cout << c << " ";
            // cout << endl;

            inner_leafs.push_back(curr_node->parent);
        }

        // cout << "I am out of loop for mean curve\n";
        // if number of leafs is odd, 
        if (num_iterations % 2 != 0) {
            curr_node = this->leafs.at(num_iterations-1); // get the last leaf

            curr_node->parent->curve = new vector<T>();
            for (auto c : *curr_node->curve) {
                curr_node->parent->curve->push_back(c); // insert the curve to the parent
            }

            inner_leafs.push_back(curr_node->parent);   // insert parent to the vector of leafs
        }



        // for (int i = 0; i < inner_leafs.size(); i++ )
        // {
        //     cout << "&&&&&&&&&&&&&&&&& " << i << " &&&&&&&&&&&&&&&&&&&&&&&\n";
        //     cout << "&&&&&&&&&&&&&&&&& SIZE: " << inner_leafs.at(i)->curve->size() << " &&&&&&&&&&&&&&&\n";
        //     for (auto c : *inner_leafs.at(i)->curve)
        //         cout << c << " ";
        //     cout << endl;
        // }


        num_iterations = inner_leafs.size();
        num_iterations_temp = num_iterations;
        this->leafs = std::move(inner_leafs);  // we go to the next level
    }


    for(auto c : *this->root->curve) {
        mean_curve.push_back(c);
    }
}


template <typename T>
void BinaryTree<T>::insert(TreeNode* node)
{
    if (this->root == NULL) {
        root = node;
        return;
    }
    else {
        std::queue<TreeNode*> queue;
        TreeNode* curr_node = this->root;

        while(curr_node != NULL) {
            if (curr_node->left == NULL) {
                curr_node->left = node;
                node->parent = curr_node;
                return;
            }
            else if (curr_node->right == NULL) {
                curr_node->right = node;
                node->parent = curr_node;
                return;
            }
            else {
                queue.push(curr_node->left);
                queue.push(curr_node->right);
            }
            curr_node = queue.front(); 
            queue.pop();
        }
    }
}

// check if a node is leaf

template <typename T>
bool BinaryTree<T>::isLeaf(TreeNode* node)
{
    if ((node->right == NULL) && (node->left == NULL))
        return true;
    return false;
}


// fill the vector of leafs
template <typename T>
void BinaryTree<T>::init_leafs(TreeNode* curr_node)
{   
    if(curr_node == NULL)
        return;

    if (isLeaf(curr_node)) {
        this->leafs.push_back(curr_node);
    }
    
    init_leafs(curr_node->left);
    init_leafs(curr_node->right);
}


// count the number of leafs which has curve
template <typename T>
int BinaryTree<T>::get_number_of_leafs()
{
    TreeNode* curr_node;
    int counter = 0;

    for (int i = 0; i < this->leafs.size(); i++)
    {
        curr_node = this->leafs.at(i);
        if (curr_node->curve != NULL)
            counter++;
    }

    return counter;
}




// this is a down to up insertion
template <typename T>
void BinaryTree<T>::insert_curve(vector<T>* curve)
{
    vector<T>* new_curve;

    for (auto leaf : this->leafs) 
    {
        if (leaf->curve == NULL) {
            // cout << "Inserting in node: " << leaf->index_node << endl;
            new_curve = new vector<T>();
            for (auto c : *curve) {
                new_curve->push_back(c);
            }
            leaf->curve = new_curve;
            return;
        }

    }
}

// Function to print binary tree in 2D
// It does reverse inorder traversal
template <typename T>
void BinaryTree<T>::print2DUtil(TreeNode *root, int space)
{
    // Base case
    if (root == NULL)
        return;
 
    // Increase distance between levels
    space += COUNT;
 
    // Process right child first
    print2DUtil(root->right, space);
 
    // Print current node after space
    // count
    cout << endl;
    for (int i = COUNT; i < space; i++)
        cout<<" ";
    cout << root->index_node <<"\n";
 
    // Process left child
    print2DUtil(root->left, space);
}


template <typename T>
void BinaryTree<T>::print_leafs()
{
    for (auto l : this->leafs)
    {
        cout << l->index_node << " ";
    }
}