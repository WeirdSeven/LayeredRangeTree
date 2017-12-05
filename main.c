#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <time.h>

typedef struct {
    double x;
    double y;
} Point;

typedef struct YArrayEntry {
    Point p;                   // the point
    struct YArrayEntry *left;  // points to y array of the left child
    struct YArrayEntry *right; // points to y array of the right child
} YArrayEntry;

typedef struct {
    double left_max;           // max of the left subtree
    double right_max;          // max of the right subtree
    YArrayEntry *y_array;      // y array
    size_t size;               // size of y array
} TreeNode;

typedef struct {
    TreeNode *nodes;
    size_t size;
} Tree;

typedef struct {
    Point *points;
    size_t size;
} PointArray;

void print_point(Point point) {
    double x = point.x;
    double y = point.y;

    if (x == DBL_MAX && y == DBL_MAX) {
        printf("(INFINITY, INFINITY)\n");
    } else if (x == DBL_MAX) {
        printf("(INFINITY, %f)\n", y);
    } else if (y == DBL_MAX) {
        printf("(%f, INFINITY)\n", x);
    } else {
        printf("(%f, %f)\n", x, y);
    }
}

void print_tree_node(TreeNode tree_node) {
    double left_max = tree_node.left_max;
    double right_max = tree_node.right_max;

    if (left_max == DBL_MAX && right_max == DBL_MAX) {
        printf("left_max = [INFINITY], right_max = [INFINITY]\n");
    } else if (left_max == DBL_MAX) {
        assert(0);
        printf("left_max = [INFINITY], right_max = [%f]\n", right_max);
    } else if (right_max == DBL_MAX) {
        printf("left_max = [%f], right_max = [INFINITY]\n", left_max);
    } else {
        printf("left_max = [%f], right_max = [%f]\n", left_max, right_max);
    }

    printf("y array size = [%zu]\n", tree_node.size);
    for (size_t i = 0; i < tree_node.size; ++i) {
        print_point(tree_node.y_array[i].p);
    }
}

void print_tree(Tree tree) {
    for (size_t i = 0; i < tree.size; ++i) {
        printf("Tree node [%zu]: \n", i);
        print_tree_node(tree.nodes[i]);
        printf("\n");
    }
}

void print_point_array(PointArray point_array) {
    for (size_t i = 0; i < point_array.size; ++i) {
        print_point(point_array.points[i]);
    }
}

double pmax(double n1, double n2) {
    if (n1 == DBL_MAX && n2 == DBL_MAX) {
        return DBL_MAX;
    } else if (n1 == DBL_MAX) {
        return n2;
    } else if (n2 == DBL_MAX) {
        return n1;
    } else {
        return n1 > n2 ? n1 : n2;
    }
}

int comparePoint(const void *a, const void *b) {
    Point *point_a = (Point *)a;
    Point *point_b = (Point *)b;
    double a_x = point_a->x;
    double b_x = point_b->x;

    if (a_x < b_x) {
        return -1;
    } else if (a_x == b_x) {
        return 0;
    } else {
        return 1;
    }
}



PointArray read_data(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error: ");
        exit(EXIT_FAILURE);
    }

    size_t num_points = 0;
    fscanf(fp, "%zu", &num_points);

    Point *points = (Point *) malloc(sizeof(Point) * num_points);

    size_t ret;
    double x, y;
    size_t count = 0;
    while ((ret = fscanf(fp, "%lf %lf", &x, &y)) != EOF) {
        points[count].x = x;
        points[count].y = y;
        count += 1;
    }

    fclose(fp);

    if (num_points != count) {
        printf("Error: Wrong number of points.\n");
        exit(EXIT_FAILURE);
    }

    qsort(points, num_points, sizeof(Point), comparePoint);

    PointArray result = { points , num_points };
    return result;
}

Tree construct_x_tree(PointArray input) {
    size_t num_points = input.size;
    size_t height = ceil(log2(num_points));
    size_t num_leaves = exp2(height);
    size_t num_nodes = num_leaves * 2;

    TreeNode *nodes = (TreeNode *) malloc(sizeof(TreeNode) * num_nodes);

    // Constructs real leave nodes
    for (size_t i = 0; i < num_points; ++i) {
        TreeNode *cur_node = &nodes[i + num_leaves];
    
        YArrayEntry *cur_y_array = (YArrayEntry *) malloc(sizeof(YArrayEntry));
        cur_y_array->p = input.points[i];
        cur_y_array->left = cur_y_array->right = NULL;

        cur_node->y_array = cur_y_array;
        cur_node->size = 1;
        cur_node->left_max = cur_node->right_max = input.points[i].x;
    }
    // Constructs dummy leave nodes
    for (size_t i = num_points; i < num_leaves; ++i) {
        TreeNode *cur_node = &nodes[i + num_leaves];

        cur_node->y_array = NULL;
        cur_node->size = 0;
        cur_node->left_max = cur_node->right_max = DBL_MAX;
    }

    // Constructs internal nodes for all levels
    for (size_t i = 1; i <= height; ++i) {
        size_t start = num_leaves >> i;
        size_t end = start << 1;

        // Constructs internal nodes with depth i
        for (size_t j = start; j < end; ++j) {
            TreeNode *cur_node = &nodes[j];
            TreeNode *left_child = &nodes[j << 1];
            TreeNode *right_child = left_child + 1;

            size_t left_size = left_child->size;
            size_t right_size = right_child->size;
            size_t cur_size = left_size + right_size;

            // Construct the y array by merging
            YArrayEntry *cur_y_array;
            // Case I: Both children are empty
            if (cur_size == 0) {
                cur_y_array = NULL;
            } else {
                cur_y_array = (YArrayEntry *) malloc(sizeof(YArrayEntry) * cur_size);

                // Case II: The right child is empty
                if (right_size == 0) {
                    for (size_t k = 0; k < left_size; ++k) {
                        YArrayEntry *cur_entry = &cur_y_array[k];
                        YArrayEntry *left_entry = &left_child->y_array[k];

                        cur_entry->left = left_entry;
                        cur_entry->right = NULL;
                        cur_entry->p = left_entry->p;
                    }
                } // Case III: Both children are non-empty
                else {
                    size_t left_ptr = 0;
                    size_t right_ptr = 0;
                    for (size_t k = 0; k < cur_size; ++k) {
                        YArrayEntry *cur_entry = &cur_y_array[k];

                        if (left_ptr == left_size) {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr - 1];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            cur_entry->p = right_entry->p;
                            right_ptr += 1;
                        } else if (right_ptr == right_size) {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr - 1];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            cur_entry->p = left_entry->p;
                            left_ptr += 1;
                        } else {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            // Sort according to the y value
                            if (left_entry->p.y <= right_entry->p.y) {
                                cur_entry->p = left_entry->p;
                                left_ptr += 1;
                            } else {
                                cur_entry->p = right_entry->p;
                                right_ptr += 1;
                            }
                        }
                    }
                }
            }

            cur_node->y_array = cur_y_array;
            cur_node->size = cur_size;
            cur_node->left_max = pmax(left_child->left_max, left_child->right_max);
            cur_node->right_max = pmax(right_child->left_max, right_child->right_max);
        }
    }

    Tree result = { nodes, num_nodes };
    return result;
}


typedef struct {
    size_t index;
    YArrayEntry *lower;
    YArrayEntry *upper;
} QueryRecord;

typedef struct {
    YArrayEntry *lower;
    YArrayEntry *upper;
    size_t size;
} QueryResultPart;

typedef struct {
    QueryResultPart *results;
    size_t size;
} QueryResult;

YArrayEntry *binary_search_lower(YArrayEntry *y_array, size_t size, double target) {
    int low = 0;
    int high = size - 1;
    YArrayEntry *found = NULL;

    while (low <= high) {
        int mid = (low + high) >> 1;

        if (y_array[mid].p.y < target) {
            low = mid + 1;
        } else {
            high = mid - 1;
            found = &y_array[mid];
        }
    }

    return found;
}

YArrayEntry *binary_search_upper(YArrayEntry *y_array, size_t size, double target) {
    int low = 0;
    int high = size - 1;
    YArrayEntry *found = NULL;

    while (low <= high) {
        int mid = (low + high) >> 1;

        if (y_array[mid].p.y > target) {
            high = mid - 1;
        } else {
            low = mid + 1;
            found = &y_array[mid];
        }
    }

    return found;
}

void follow_y_array_entry_pointers(YArrayEntry **lower_ptr, YArrayEntry **upper_ptr, 
                                    double lower_value, double upper_value,
                                    YArrayEntry *lower_pte, YArrayEntry*upper_pte) {
    if (lower_value > lower_pte->p.y) {
        *lower_ptr = lower_pte + 1;
    } else {
        *lower_ptr = lower_pte;
    }

    if (upper_value < upper_pte->p.y) {
        *upper_ptr = upper_pte - 1;
    } else {
        *upper_ptr = upper_pte;
    }
}

QueryResult query(Tree tree, Point lower_left, Point upper_right) {
    struct timespec begin, end, end2;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    double x_lower = lower_left.x;
    double x_upper = upper_right.x;
    double y_lower = lower_left.y;
    double y_upper = upper_right.y;

    QueryResult null_result = { NULL, 0};
    if (x_lower > x_upper) {
        printf("x lower is greater than x upper!\n");
        return null_result;
    }
    if (y_lower > y_upper) {
        printf("y lower is greater than y upper!\n");
        return null_result;
    }

     // This is the textbook definition of height plus 1
    size_t height = log2(tree.size);

    QueryRecord *x_lower_path_record = (QueryRecord *) malloc(sizeof(QueryRecord) * height);
    QueryRecord *x_upper_path_record = (QueryRecord *) malloc(sizeof(QueryRecord) * height);

    TreeNode *root = &tree.nodes[1];
    YArrayEntry *y_lower_ptr_root = binary_search_lower(root->y_array, root->size, y_lower);
     //printf("Here\n");
    if (y_lower_ptr_root == NULL) {
        printf("y lower is greater than max!\n");
        return null_result;
    }
    YArrayEntry *y_upper_ptr_root = binary_search_upper(root->y_array, root->size, y_upper);
    if (y_upper_ptr_root == NULL) {
        printf("y upper is less than min!\n");
        return null_result;
    }

     // Search for x_lower
    size_t cur_node_index = 1;
    YArrayEntry *y_lower_ptr = y_lower_ptr_root;
    YArrayEntry *y_upper_ptr = y_upper_ptr_root;
    int x_lower_gt_max = 0;
    size_t x_lower_search_height;

    for (size_t i = 0; i < height; ++i) {
        x_lower_search_height = i;
        QueryRecord *cur_record = &x_lower_path_record[i];
        TreeNode *cur_node = &tree.nodes[cur_node_index];

        // We know if x_lower is greater than max
        // exactly at the root level
        if (x_lower > cur_node->right_max) {
            x_lower_gt_max = 1;
            break;
        }

        // The leave level is handled differently
        if (i == height - 1) {
            assert(x_lower <= cur_node->left_max);
            assert(y_lower_ptr == y_upper_ptr);
            cur_record->index = cur_node_index;
            cur_record->lower = y_lower_ptr;
            cur_record->upper = y_upper_ptr;      
        } else {
            double y_lower_value = y_lower_ptr->p.y;
            double y_upper_value = y_upper_ptr->p.y;
            YArrayEntry *y_lower_ptr_right = y_lower_ptr->right;
            YArrayEntry *y_upper_ptr_right = y_upper_ptr->right;
            YArrayEntry *y_lower_ptr_left = y_lower_ptr->left;
            YArrayEntry *y_upper_ptr_left = y_upper_ptr->left;

            // Take a left
            if (x_lower <= cur_node->left_max) {
                cur_record->index = cur_node_index;
                follow_y_array_entry_pointers(&cur_record->lower, &cur_record->upper,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_right, y_upper_ptr_right);
                // There are no corresponding points in this subtree
                if (cur_record->lower > cur_record->upper) {
                    cur_record->lower = cur_record->upper = NULL;
                }


                cur_node_index = cur_node_index << 1;
                follow_y_array_entry_pointers(&y_lower_ptr, &y_upper_ptr,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_left, y_upper_ptr_left);
                // There are no corresponding points in this subtree. Stop searching
                if (y_lower_ptr > y_upper_ptr) {
                    break;
                }
            } // Take a right
            else {
                cur_record->index = cur_node_index;
                cur_record->lower = NULL;
                cur_record->upper = NULL;

                cur_node_index = (cur_node_index << 1) + 1;
                follow_y_array_entry_pointers(&y_lower_ptr, &y_upper_ptr,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_right, y_upper_ptr_right);
                if (y_lower_ptr > y_upper_ptr) {
                    break;
                }
            } 
        }
    }

    if (x_lower_gt_max) {
        printf("x lower is greater than max!\n");
        return null_result;
    }

     // Search for x_upper
    cur_node_index = 1;
    y_lower_ptr = y_lower_ptr_root;
    y_upper_ptr = y_upper_ptr_root;
    int x_upper_lt_min = 0;
    size_t x_upper_search_height;

    for (size_t i = 0; i < height; ++i) {
        x_upper_search_height = i;
        QueryRecord *cur_record = &x_upper_path_record[i];
        TreeNode *cur_node = &tree.nodes[cur_node_index];

        // We only know if x_upper is less than min
        // at the leave level. Here left_max = right_max
        if (cur_node_index == tree.size >> 1 && x_upper < cur_node->left_max) {
            x_upper_lt_min = 1;
            break;
        }

        if (i == height - 1) {
            assert(y_lower_ptr == y_upper_ptr);
            cur_record->index = cur_node_index;
            if (x_upper >= cur_node->left_max) {
                cur_record->lower = y_lower_ptr;
                cur_record->upper = y_upper_ptr; 
            } else {
                cur_record->lower = cur_record->upper = NULL;
            } 
        } else {
            double y_lower_value = y_lower_ptr->p.y;
            double y_upper_value = y_upper_ptr->p.y;
            YArrayEntry *y_lower_ptr_right = y_lower_ptr->right;
            YArrayEntry *y_upper_ptr_right = y_upper_ptr->right;
            YArrayEntry *y_lower_ptr_left = y_lower_ptr->left;
            YArrayEntry *y_upper_ptr_left = y_upper_ptr->left;

            // Take a left. Be careful to handle duplicates correctly
            // and not walk off the tree
            if ((x_upper <= cur_node->left_max && x_upper < cur_node->right_max)
                || cur_node->right_max == DBL_MAX) {
                //printf("Take a left\n");
                cur_record->index = cur_node_index;
                cur_record->lower = NULL;
                cur_record->upper = NULL;

                cur_node_index = cur_node_index << 1;
                follow_y_array_entry_pointers(&y_lower_ptr, &y_upper_ptr,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_left, y_upper_ptr_left);

                if (y_lower_ptr > y_upper_ptr) {
                    break;
                }
            } // Take a right
            else {
                cur_record->index = cur_node_index;
                follow_y_array_entry_pointers(&cur_record->lower, &cur_record->upper,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_left, y_upper_ptr_left);
                if (cur_record->lower > cur_record->upper) {
                    cur_record->lower = cur_record->upper = NULL;
                }

                cur_node_index = (cur_node_index << 1) + 1;
                follow_y_array_entry_pointers(&y_lower_ptr, &y_upper_ptr,
                                                y_lower_value, y_upper_value,
                                                y_lower_ptr_right, y_upper_ptr_right);
                if (y_lower_ptr > y_upper_ptr) {
                    break;
                }
            }
        }
    }

    if (x_upper_lt_min) {
        printf("x upper is less than min!\n");
        return null_result;
    }



    // Search for the first level where the two paths differ
    size_t split_height = height;
    {
        size_t i = 0;
        while (i <= x_lower_search_height && i <= x_upper_search_height) {
            if (x_lower_path_record[i].index == x_upper_path_record[i].index) {
                ++i;
            } else {
                break;
            }
        }
        split_height = i;
    }

    QueryResultPart *x_lower_path_result = (QueryResultPart *) malloc(sizeof(QueryResultPart) * (height - split_height));
    QueryResultPart *x_upper_path_result = (QueryResultPart *) malloc(sizeof(QueryResultPart) * (height - split_height));

    // Collect results for the lower branch
    size_t x_lower_result_size = 0;
    for (size_t i = x_lower_search_height; i >= split_height; --i) {
        QueryRecord *cur_record = &x_lower_path_record[i];

        if (cur_record->lower && cur_record->upper) {
            QueryResultPart *cur_result = &x_lower_path_result[x_lower_result_size];
            cur_result->lower = cur_record->lower;
            cur_result->upper = cur_record->upper;
            cur_result->size = cur_record->upper - cur_record->lower + 1;

            x_lower_result_size += 1;
        }
    }

    // Collect results for the upper branch
    size_t x_upper_result_size = 0;
    for (size_t i = split_height; i <= x_upper_search_height; ++i) {
        QueryRecord *cur_record = &x_upper_path_record[i];

        if (cur_record->lower && cur_record->upper) {
            QueryResultPart *cur_result = &x_upper_path_result[x_upper_result_size];
            cur_result->lower = cur_record->lower;
            cur_result->upper = cur_record->upper;
            cur_result->size = cur_record->upper - cur_record->lower + 1;

            x_upper_result_size += 1;
        }
    }

    size_t result_size = x_lower_result_size + x_upper_result_size ;
    QueryResultPart *results = (QueryResultPart *) malloc(sizeof(QueryResultPart) * result_size);
    memcpy(results, x_lower_path_result, sizeof(QueryResultPart) * x_lower_result_size);
    memcpy(results + x_lower_result_size, x_upper_path_result, sizeof(QueryResultPart) * x_upper_result_size);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec  - begin.tv_sec);
    printf("Range Tree Time: %f\n", elapsed_time);

    QueryResult result = { results, result_size };
    return result;
}

size_t get_num_points_from_query_result(QueryResult query_result) {
    size_t size = 0;
    for (size_t i = 0; i < query_result.size; ++i) {
        size += query_result.results[i].size;
    }
    return size;
}

PointArray get_point_array_from_query_result(QueryResult query_result) {
    size_t num_points = get_num_points_from_query_result(query_result);

    Point *points = (Point *) malloc(sizeof(Point) * num_points);
    size_t cur_index = 0;
    for (size_t i = 0; i < query_result.size; ++i) {
        QueryResultPart *cur_result = &query_result.results[i];
        for (size_t j = 0; j < cur_result->size; ++j) {
            YArrayEntry *cur_entry = cur_result->lower + j;
            points[cur_index] = cur_entry->p;
            cur_index += 1;
        }
    }

    PointArray result = { points, num_points };
    return result;
}

PointArray linear_query(PointArray point_array, Point lower_left, Point upper_right) {
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    Point *points = (Point *) malloc(sizeof(Point) * point_array.size);
    size_t num_points = 0;
    for (size_t i = 0; i < point_array.size; i++) {
        Point cur_point = point_array.points[i];
        if (cur_point.x >= lower_left.x && cur_point.y >= lower_left.y 
            && cur_point.x <= upper_right.x && cur_point.y <= upper_right.y) {
            points[num_points] = cur_point;
            num_points += 1;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed_time = (end.tv_nsec - begin.tv_nsec) / 1000000000.0 + (end.tv_sec  - begin.tv_sec);
    printf("Linear Time: %f\n", elapsed_time);

    PointArray result = { points, num_points };
    return result;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Usage: %s datafile\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    PointArray input = read_data(argv[1]);

    Tree x_tree = construct_x_tree(input);
    printf("Done building the tree!\n");

    //print_tree(x_tree);

    do {
        Point lower_left;
        Point upper_right;

        printf("Please enter the query (in the order of lower-left.x lower-left.y upper-right.x upper-right.y, enter all zeros to terminate): ");
        scanf("%lf %lf %lf %lf", &lower_left.x, &lower_left.y, &upper_right.x, &upper_right.y);

        if (lower_left.x == 0 && lower_left.y == 0 & upper_right.x == 0 && upper_right.y == 0)
            break;

        QueryResult query_result = query(x_tree, lower_left, upper_right);
        PointArray point_array = get_point_array_from_query_result(query_result);
        printf("Range: \n");
        //print_point_array(point_array);
        PointArray linear_point_array = linear_query(input, lower_left, upper_right);
        printf("Linear: \n");
        //print_point_array(linear_point_array);

        free(point_array.points);
        free(linear_point_array.points);
    } while (1);

    return 0;
}