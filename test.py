import unittest

from tree import *


class TestTree(unittest.TestCase):
    
    def test_all_trees(self):
        self.assertEqual(set(all_trees([1,2])),
            set([
                (1,2),
        ]))
        
        self.assertEqual(set(all_trees([1,3,7])),
            set([
                (1,3,7),
                (1,(3,7)),
                (3,(1,7)),
                (7,(1,3)),
        ]))
        
        self.assertEqual(set(all_trees([1,2,3,4])),
            set([
                (1,2,3,4),
                (1,2,(3,4)),
                (1,3,(2,4)),
                (1,4,(2,3)),
                (2,3,(1,4)),
                (2,4,(1,3)),
                (3,4,(1,2)),
                (1,(2,3,4)),
                (2,(1,3,4)),
                (3,(1,2,4)),
                (4,(1,2,3)),
                (1,(2,(3,4))),
                (1,(3,(2,4))),
                (1,(4,(2,3))),
                (2,(1,(3,4))),
                (2,(3,(1,4))),
                (2,(4,(1,3))),
                (3,(1,(2,4))),
                (3,(2,(1,4))),
                (3,(4,(1,2))),
                (4,(1,(2,3))),
                (4,(2,(1,3))),
                (4,(3,(1,2))),
                ((1,2),(3,4)),
                ((1,3),(2,4)),
                ((1,4),(2,3)),
        ]))
    
    
    def test_leaf_set(self):
        self.assertEqual(leaf_set((1,2)), [1,2])
        self.assertEqual(leaf_set((1,2,3)), [1,2,3])
        self.assertEqual(leaf_set((1,(2,3))), [1,2,3])
        self.assertEqual(leaf_set((2,5,9)), [2,5,9])
        self.assertEqual(leaf_set((1,(3,4,(2,6)))), [1,2,3,4,6])
        self.assertEqual(leaf_set(((1,4),(2,5))), [1,2,4,5])
    
    
    def test_clusters(self):
        for i in xrange(2):
            self.assertEqual(set(clusters((1,2))), set([
                (1,),
                (2,),
                (1,2),
            ]))
            
            self.assertEqual(set(clusters((1,(2,3)))), set([
                (1,),
                (2,),
                (3,),
                (2,3),
                (1,2,3),
            ]))
            
            self.assertEqual(set(clusters((1,2,3,4))), set([
                (1,),
                (2,),
                (3,),
                (4,),
                (1,2,3,4),
            ]))
            
            self.assertEqual(set(clusters((1,2,(3,4)))), set([
                (1,),
                (2,),
                (3,),
                (4,),
                (3,4),
                (1,2,3,4),
            ]))
            
            self.assertEqual(set(clusters((1,((2,3),(4,5))))), set([
                (1,),
                (2,),
                (3,),
                (4,),
                (5,),
                (2,3),
                (4,5),
                (2,3,4,5),
                (1,2,3,4,5),
            ]))
    
    
    def test_compare(self):
        for i in xrange(2):
            self.assertTrue(compare((1,2,3), (1,(2,3))))
            self.assertTrue(compare((1,2,3), (2,(1,3))))
            self.assertTrue(compare((1,2,3), (3,(1,2))))
            self.assertTrue(compare((1,2,3), (1,2,3)))
            self.assertFalse(compare((1,(2,3)), (1,2,3)))
            self.assertFalse(compare((2,(1,3)), (1,2,3)))
            self.assertFalse(compare((3,(1,2)), (1,2,3)))
            
            self.assertTrue(compare((1,(2,3,(4,5))), (1,((2,3),(4,5)))))
            self.assertTrue(compare((1,(2,3,(4,5))), (1,(2,(3,(4,5))))))
            self.assertFalse(compare((1,(2,3,(4,5))), (1,(2,(3,4,5)))))
    
    
    def test_restriction(self):
        self.assertEqual(restriction((1,2,3), [1,4,5]), 1)
        self.assertEqual(restriction((1,2,3), [4,5]), None)
        self.assertEqual(restriction((1,2,3,4), [2,3,4]), (2,3,4))
        self.assertEqual(restriction((1,(2,3,4)), [2,3,4]), (2,3,4))
        self.assertEqual(restriction((1,(2,3,4)), [1,2,3]), (1,(2,3)))
        self.assertEqual(restriction((1,(2,3,(4,5))), [1,2,4,5]), (1,(2,(4,5))))
        self.assertEqual(restriction((1,(2,3,(4,5))), [1,2,3,4,5]), (1,(2,3,(4,5))))


if __name__ == '__main__':
    unittest.main()
