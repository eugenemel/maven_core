/****************************************************************************

    Neural Network Library
    Copyright (C) 1998 Daniel Franklin

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
    02111-1307, USA.
              
 ****************************************************************************/

// This library implements a very basic three layer backpropagation neural
// network. The two classes are "neuron" which represents an adaptive line
// combiner similar to ADALINE, and "nnwork" which is the entire network.

#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <neuron.h>


// Constructor.

nnlayer::nnlayer (int dimension, int num_weights)
{
	size = dimension;
	weights = num_weights;
    nodes = std::vector<neuron>(static_cast<unsigned long>(size)); //OLD: new neuron [size];
//OLD: assert (nodes);
	
	for (int j = 0; j < size; j++) {
        nodes.at(j).weights = std::vector<float>(static_cast<unsigned long>(weights));//OLD nodes [j].weights = new float [weights];
        //OLD assert (nodes [j].weights);

        for (unsigned int i = 0; i < weights; i++)
            nodes.at(j).weights.at(i) = 0.5 - rand() / float (RAND_MAX);
	}
}

//OLD: destructor necessary to clean up pointers.
// Destructor. 

//nnlayer::~nnlayer ()
//{
//	for (int i = 0; i < size; i++)
//		delete [] nodes [i].weights;

//	delete [] nodes;
//}
