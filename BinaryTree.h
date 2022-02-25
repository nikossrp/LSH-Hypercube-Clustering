#ifndef __BinaryTree__
#define __BinaryTree__
#include "Utilities.h"
// #include "Utilities.cpp"
#define COUNT 10


typedef struct TreeNode {
    std::vector<Type>* curve;         // The data in this node.
    int index_node;
    struct TreeNode *left;   // Pointer to the left subtree.
    struct TreeNode *right;  // Pointer to the right subtree.
    struct TreeNode *parent; 

}TreeNode;


template <typename T>
class BinaryTree {
    private:
        TreeNode* root;
        int num_leafs;
        int num_nodes;
        int height;
        std::vector<TreeNode*> leafs;       // we need that for inserting the curves

        void insert(TreeNode* node);        // insert at the first available leaf
        void init_leafs(TreeNode* root);    // fill the vector of leafs
        void destroy_tree(TreeNode* root);
        bool isLeaf(TreeNode* node);

    public:
        BinaryTree(int num_curves);
        ~BinaryTree();

        void insert_curve(std::vector<T>* curve);
        void get_mean_curve(std::vector<T>& mean_curve);    // compute the mean curve, - mean curve will be in the root at the end of caclculations -
        int get_number_of_leafs();      // count the number of leafs which has curve

        // for debugging
        TreeNode* get_root() { return root; }
        void print2DUtil(TreeNode *root, int space);
        void print_leafs();

};

#endif 

