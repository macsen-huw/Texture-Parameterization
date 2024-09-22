///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//  
///////////////////////////////////////////////////

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <queue>
#include <cmath>
#include <map>
#include <algorithm>
// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))

#define N_ITERATIONS 100000

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
	otherHalf.resize(0);

    numBoundary = 0;
    } // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;
                
            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();
                
                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        geometryStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;                  
                        } // tex coord
                    default:
                        break;
                    } // switch on second character 
                break;
                } // some sort of vertex data
                
            case 'f':       // face data
                { // face
				// make a hard assumption that we have a single triangle per line
                unsigned int vertexID;
                
                // read in three vertices
				for (unsigned int vertex = 0; vertex < 3; vertex++)
					{ // per vertex
					// read a vertex ID
					geometryStream >> vertexID;

					// subtract one and store them (OBJ uses 1-based numbering)
					faceVertices.push_back(vertexID-1);
					} // per vertex
				break;
                } // face
                
            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
                
            } // per vertex
        } // non-empty vertex set

// 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
// 	std::cout << "Object Size:       " << objectSize << std::endl;

    // return a success code
    return true;
	} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex].x << " " << textureCoords[vertex].y << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        geometryStream << "f";
        
        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
			} // per vertex
		// end the line
        geometryStream << std::endl;
        } // per face
    
    } // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
    { // Render()
	// make sure that textures are disabled
	glDisable(GL_TEXTURE_2D);

	float scale = renderParameters->zoomScale;
	scale /= objectSize;
	// Scale defaults to the zoom setting

    //Only move the image if showing the mesh, not the texture
    if(!(renderParameters->renderTexture || renderParameters->renderNormalMap))
        glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);
		
	if (renderParameters->useWireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // start rendering
    glBegin(GL_TRIANGLES);
	
    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        
		// now do a loop over three vertices
		for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex

            if(renderParameters->renderNormalMap)
            {
                glColor3f(normals[faceVertices[face+vertex]].x, normals[faceVertices[face+vertex]].y, normals[faceVertices[face+vertex]].z);

                glVertex3f
                    (
                    (textureCoords[faceVertices[face+vertex]].x - 0.5) * 2,
                    (textureCoords[faceVertices[face+vertex]].y - 0.5) * 2,
                    textureCoords[faceVertices[face+vertex]].z
                    );

            }

            else
            {
                if(renderParameters->useNormal)
                    glColor3f(normals[faceVertices[face+vertex]].x, normals[faceVertices[face+vertex]].y, normals[faceVertices[face+vertex]].z);

                else if(renderParameters->useTexCoords)
                    glColor3f(textureCoords[faceVertices[face+vertex]].x, textureCoords[faceVertices[face+vertex]].y, textureCoords[faceVertices[face+vertex]].z);

                else{
                    // set colour using vertex ID
                    glColor3f
                        (
                        colours[faceVertices[face+vertex]].x,
                        colours[faceVertices[face+vertex]].y,
                        colours[faceVertices[face+vertex]].z
                        );

                }

                if(renderParameters->renderTexture)
                {
                    glVertex3f
                        (
                        (textureCoords[faceVertices[face+vertex]].x - 0.5) * 2,
                        (textureCoords[faceVertices[face+vertex]].y - 0.5) * 2,
                        textureCoords[faceVertices[face+vertex]].z
                        );
                }

                else
                {
                    // use scaled xyz for vertex position
                    glVertex3f
                        (
                        scale * vertices[faceVertices[face+vertex]].x,
                        scale * vertices[faceVertices[face+vertex]].y,
                        scale * vertices[faceVertices[face+vertex]].z
                        );
                }

            }

			} // per vertex
        } // per face

    // close off the triangles
    glEnd();

    // revert render mode  
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    } // Render()

void AttributedObject::setup()
{
    //We already have vertex and face data
    //For a directed edge structure, we also need the first directed edge and the other half
    //First of all, we initialise firstDirectedEdge vector (for simplicity)
    for (int i = 0; i < vertices.size(); i++)
    {
        firstDirectedEdge.push_back(-1);
    }

    //Initialise other half vector
    for (int i = 0; i < faceVertices.size(); i++)
    {
        otherHalf.push_back(-1);
    }

    //Loop through all vertices in the face
    for (int i = 0; i < faceVertices.size(); i++)
    {

        //Has the first directed edge been placed?
        if (firstDirectedEdge[faceVertices[i]] == -1)
        {
            //Find the edge that goes towards it (since edges are pointing "to")
            int faceNumber = i / 3;
            int faceIndex = i % 3;
            firstDirectedEdge[faceVertices[i]] = 3*faceNumber + ((faceIndex + 1) % 3);
        }

        int edgeFace = i / 3;
        int edgeFaceIndex = i % 3;
        int vertexFrom = faceVertices[3*edgeFace + ((edgeFaceIndex + 2) % 3)];
        int vertexTo = faceVertices[3*edgeFace + (edgeFaceIndex % 3)];

        //Loop through all faces again to find other half
        for(int j = 0; j < faceVertices.size(); j+= 3)
        {


            if(faceVertices[j] == vertexTo && faceVertices[j+1] == vertexFrom)
            {
                otherHalf[i] = j+1;
                break;
            }

            if(faceVertices[j+1] == vertexTo && faceVertices[j+2] == vertexFrom)
            {
                otherHalf[i] = j+2;
                break;
            }

            if(faceVertices[j+2] == vertexTo && faceVertices[j] == vertexFrom)
            {
                otherHalf[i] = j;
                break;
            }
        }
    }


    //Now that we've set up the data structure, we perform Floater's algorithm
    performFloater();

    //After completing this, we calculate vertex normals
    if(normals.empty())
            calculateNormals();

    //Get the max and min values of different parts
    float maxTexX, maxTexY, minTexX, minTexY;
    for(Cartesian3 tex : textureCoords)
    {
        if(tex.x > maxTexX)
            maxTexX = tex.x;
        if(tex.y > maxTexY)
            maxTexY = tex.y;
        if(tex.x < minTexX)
            minTexX = tex.x;
        if(tex.y < minTexY)
            minTexY = tex.y;
    }

    //std::cout << "Highest values: " << maxTexX << ", " << maxTexY << std::endl;
    //std::cout << "Lowest values: " << minTexX << ", " << minTexY << std::endl;

}

//Find the boundary vertices
void AttributedObject::findBoundary()
{

    //Given the truncation, we only want the vertices that are used in the shape (i.e. we don't want the vertices that show Hamish's neck)
    //We find these by seeing which vertices are within the faces
    std::vector<int> usedVertices(vertices.size(), 0);

   //Find all vertices being used
    for(int vertex: faceVertices)
    {
        if(usedVertices[vertex] == 0)
            usedVertices[vertex] = 1;
    }


    //Finding the boundary vertices of the patch
    //In other words, all edges without other halves
    for (int i = 0; i < otherHalf.size(); i++)
    {
        if(otherHalf[i] == -1)
        {
            //This edge is on the boundary,
            //Find out the vertex that it has come from
            int edgeFace = i / 3;
            int originalIndex = (i % 3) + 2;

            int boundaryVertex = faceVertices[3*edgeFace + (originalIndex % 3)];

            //Find out where the edge was going to
            int nextBoundaryVertex = faceVertices[3*edgeFace + ((originalIndex + 1) % 3)];

            //Test whether the boundaryVertex is being used
            if(usedVertices[boundaryVertex] == 1)
            {
                //Add to the vertices if so
                boundaryVertices[boundaryVertex] = 1;
                boundaryConnections[boundaryVertex] = (nextBoundaryVertex);

                //Remove the vertex from the usedVertex (makes the next step easier, also removes duplicates)
                usedVertices[boundaryVertex] = 0;
            }

        }

    }

    //Finally, get the interior vertices (i.e. any used vertices that aren't boundaries)
    for (int i = 0; i < usedVertices.size(); i++)
    {
        if(usedVertices[i] == 1)
            interiorVertices[i] = 1;
    }
}

void AttributedObject::createPatches()
{

    //Learn how many vertices are being used in the first place
    std::vector<int> usedVertices(vertices.size(), 0);

    for (int i = 0; i < firstDirectedEdge.size(); i++)
    {
        if(firstDirectedEdge[i] != -1)
            usedVertices[i] = 1;
    }

    patchIndex.resize(vertices.size(), -1);

    //Initialise important variables
    std::deque<int> queue;
    int patchSize = 0; //Determines how many vertices are in the patch
    int patchNumber = 0; //Assigns a patch index to a vertex

    //We keep iterating until all vertices have been assigned to a patch
    bool minusNotPresent;
    do
    {
        minusNotPresent = true;

        //Empty the queue
        queue.clear();

        //Find vertex to start generation with
        for(auto pair : interiorVertices)
        {
            if(patchIndex[pair.first] == -1)
            {
                queue.push_back(pair.first);

                //Add this first element to the patch
                patchIndex[pair.first] = patchNumber;
                patchSize++;
                break;
            }

        }

        //If there isn't an interior vertex, then we select a boundary vertex
        if(queue.empty())
        {
            for(auto pair : boundaryVertices)
            {
                if(patchIndex[pair.first] == -1)
                {
                    queue.push_back(pair.first);

                    //Add this first element to the patch
                    patchIndex[pair.first] = patchNumber;
                    patchSize++;
                    break;
                }

            }

        }

        //Keep looping until queue runs out, or the patch exceeds 2000 vertices
        while(!(queue.empty()) && patchSize <= 2000)
        {

            //Take first value of queue
            int queueFront = queue.at(0);
            queue.pop_front();

            //Walk around the vertex at the front of the queue
            int startEdge = firstDirectedEdge[queueFront];
            int nextEdge = startEdge;

            do
            {

                int currentEdgeFace = nextEdge / 3;
                int currentEdgeFaceIndex = nextEdge % 3;

                int neighbour = faceVertices[3*currentEdgeFace + currentEdgeFaceIndex];

                //Add to the patch if it isn't already
                if(patchIndex[neighbour] == -1)
                {
                    patchIndex[neighbour] = patchNumber;
                    patchSize++;

                    //Check if it's already in the queue
                    bool inQueue = false;
                    for(int item : queue)
                    {
                        if(item == neighbour)
                            inQueue = true;
                    }

                    //If not, then add it
                    if (!inQueue)
                        queue.push_back(neighbour);
                }

                int otherHalfEdge = otherHalf[nextEdge];

                //Check that it isn't -1
                if(otherHalfEdge == -1)
                    break;

                //Discover the vertex belonging to the other half
                int otherEdgeFace = otherHalfEdge / 3;
                int otherEdgeIndex = otherHalfEdge % 3;

                //Go to next edge
                nextEdge = 3*otherEdgeFace + ((otherEdgeIndex + 1) % 3);

            } while (nextEdge != startEdge);

            //We've successfully walked around the vertex, so we go to the next item in the queue
        }

        //Here, we've either emptied the queue, or reached the maximum number of vertices in a patch
        //Final check - is every vertex assigned a patch?

        for(int i = 0; i < patchIndex.size(); i++)
        {
            if(patchIndex[i] == -1 && usedVertices[i] == 1)
                minusNotPresent = false;
        }

        //Iterate patch index, and reset number of vertices in a patch
        patchNumber++;
        patchSize = 0;


      //Keep repeating until every vertex is assigned a patch
    } while (minusNotPresent == false);

}

//Map the boundary edges to the outside of a square
void AttributedObject::mapBoundaryVertices()
{
    //Due to the orientation of the model, we use the y and z values to determine mapping
    float highestZ = -5;
    float lowestZ = 5;
    float highestY = -5;
    float lowestY = 5;
    float highestX = -5;
    float lowestX = 5;

    //Loop through boundary vertices to get highest and lowest values of each coordinate
    for (auto pair : boundaryVertices)
    {
        Cartesian3 vertex = vertices[pair.first];

        //Get the highest and lowest coordinate of the vertices
        highestZ = std::max(highestZ, vertex.z);
        highestY = std::max(highestY, vertex.y);
        highestX = std::max(highestX, vertex.x);
        lowestZ = std::min(lowestZ, vertex.z);
        lowestY = std::min(lowestY, vertex.y);
        lowestX = std::min(lowestX, vertex.x);
     }

    //Get the dimensions of the square that will be made

    //axesOrientation = ZY_COORDINATES;
    float squareWidth, squareHeight;

    if(axesOrientation == ZY_COORDINATES)
    {
        squareWidth = highestZ - lowestZ;
        squareHeight = highestY - lowestY;
    }

    else if (axesOrientation == XZ_COORDINATES)
    {
        squareWidth = highestX - lowestX;
        squareHeight = highestZ - lowestZ;
    }

    else
    {
        squareWidth = highestX - lowestX;
        squareHeight = highestY - lowestY;
    }


    std::vector<Cartesian3> textureCoordinates(vertices.size(), {0,0,0});

    //Get every boundary vertex
    for (auto pair : boundaryVertices)
    {
        Cartesian3 vertex = vertices[pair.first];

        //Find out what the closest distance is (i.e. nearest to 0)
        float distAcross, distUpwards;


        /*
        distAcross = std::abs(vertex.z - lowestZ);

        if(ZYCoordinates)
            distUpwards = std::abs(vertex.y - lowestY);
        else
            distUpwards = std::abs(vertex.x - lowestX);
        */

        switch (axesOrientation)
        {
            case ZY_COORDINATES:
                distAcross = std::abs(vertex.z - lowestZ);
                distUpwards = std::abs(vertex.y - lowestY);
            break;

            case XZ_COORDINATES:
                distAcross = std::abs(vertex.x - lowestX);
                distUpwards = std::abs(vertex.z - lowestZ);
            break;

            //default - XY
            default:
                distAcross = std::abs(vertex.x - lowestX);
                distUpwards = std::abs(vertex.y - lowestY);
            break;
        }



        float acrossCoordinate = distAcross / squareWidth;
        float upCoordinate = distUpwards / squareHeight;

        float closerToRight = 1 - acrossCoordinate;
        float closerToUp = 1 - upCoordinate;

        //Get the closest distance
        float closestDistance = std::min(acrossCoordinate, std::min(upCoordinate, std::min(closerToRight, closerToUp)));

        //Closest to top
        if(closestDistance == closerToUp)
            textureCoordinates[pair.first] = {acrossCoordinate, 1, 0};

        //Closest to right hand side
        else if (closestDistance == closerToRight)
            textureCoordinates[pair.first] = {1, upCoordinate, 0};

        //Closest to bottom
        else if(closestDistance == upCoordinate)
            textureCoordinates[pair.first] = {acrossCoordinate, 0, 0};

        //Closest to left hand side
        else if (closestDistance == acrossCoordinate)
            textureCoordinates[pair.first] = {0, upCoordinate, 0};

    }

    //We've now found the boundary coordinates, so we can count the number of them
    numBoundary = boundaryVertices.size();

    //We keep track of the closest boundaryCoordinate to the corner so we can clamp to the corner of the square
    //The first index is the distance, the second is the vertexID
    float topLeft[2] = {100, -1};
    float topRight[2] = {100, -1};
    float bottomLeft[2] = {100, -1};
    float bottomRight[2] = {100, -1};


    //Now loop through the boundary vertices to see which vertices belong to the corners
    for (auto pair : boundaryVertices)
    {
        Cartesian3 point = textureCoordinates[pair.first];

        //Get the distances from every corner
        float distanceToTopLeft = (Cartesian3(0,1,0) - point).length();
        float distanceToTopRight = (Cartesian3(1,1,0) - point).length();
        float distanceToBottomLeft = (Cartesian3(0,0,0) - point).length();
        float distanceToBottomRight = (Cartesian3(1,0,0) - point).length();

        //Start with top left corner
        if (distanceToTopLeft < topLeft[0])
        {
            topLeft[0] = distanceToTopLeft;
            topLeft[1] = pair.first;
        }

        //Top Right Corner
        if (distanceToTopRight < topRight[0])
        {
            topRight[0] = distanceToTopRight;
            topRight[1] = pair.first;
        }

        //Bottom Left Corner
        if (distanceToBottomLeft < bottomLeft[0])
        {
            bottomLeft[0] = distanceToBottomLeft;
            bottomLeft[1] = pair.first;
        }

        //Bottom Right Corner
        if(distanceToBottomRight < bottomRight[0])
        {
            bottomRight[0] = distanceToBottomRight;
            bottomRight[1] = pair.first;
        }

    }

    //We now know the corner variables, so we can update the coordinates
    textureCoordinates[topLeft[1]] = {0,1,0};
    textureCoordinates[topRight[1]] = {1,1,0};
    textureCoordinates[bottomLeft[1]] = {0,0,0};
    textureCoordinates[bottomRight[1]] = {1,0,0};

    //We've now fixed to the corners, so we make sure that all vertices between corners are clamped to the correct wall
    //Take the vertexIDs of the corner variables
    int topLeftVertex = topLeft[1];
    int topRightVertex = topRight[1];
    int bottomLeftVertex = bottomLeft[1];
    int bottomRightVertex = bottomRight[1];


    //Go from bottom right to top right, ensuring that all values are correctly mapped to the top
    int next = bottomRightVertex;
    do
    {
        next = boundaryConnections[next];

        //Check whether the boundary coordinate is clamped to the right hand edge
        if (textureCoordinates[next][0] != 1)
        {
            //Redo the barycentric mapping, ensuring it's to the correct edge
            Cartesian3 vertex = vertices[next];

            float distToDown;
            if(axesOrientation == XZ_COORDINATES)
                distToDown = std::abs(vertex.z - lowestZ);
            else
                distToDown = std::abs(vertex.y - lowestY);

            float yCoordinate = distToDown / squareHeight;

            //Update the coordinate
            textureCoordinates[next] = {1, yCoordinate, 0};

        }

    } while(next != topRightVertex);

    //Top Right -> Top Left
    next = topRightVertex;
    do
    {
        next = boundaryConnections[next];

        //Check whether the boundary coordinate is clamped to the upper edge
        if (textureCoordinates[next][1] != 1)
        {
            //Redo the barycentric mapping, ensuring it's to the correct edge
            Cartesian3 vertex = vertices[next];

            float distToLeft;
            if(axesOrientation == ZY_COORDINATES)
                distToLeft = std::abs(vertex.z - lowestZ);
            else
                distToLeft = std::abs(vertex.x - lowestX);

            float zCoordinate = distToLeft / squareWidth;

            //Update the coordinate
            textureCoordinates[next] = {zCoordinate, 1, 0};
        }

    } while(next != topLeftVertex);

    //Top Left -> Bottom Left
    next = topLeftVertex;
    do
    {
        next = boundaryConnections[next];

        //Check whether the boundary coordinate is clamped to the left hand edge
        if (textureCoordinates[next][0] != 0)
        {
            //Redo the barycentric mapping, ensuring it's to the correct edge
            Cartesian3 vertex = vertices[next];

            float distToDown;

            if(axesOrientation == XZ_COORDINATES)
                distToDown = std::abs(vertex.z - lowestZ);
            else
                distToDown = std::abs(vertex.y - lowestY);

            float yCoordinate = distToDown / squareHeight;

            //Update the coordinate
            textureCoordinates[next] = {0, yCoordinate, 0};

        }

    } while(next != bottomLeftVertex);

    //Bottom Left -> Bottom Right
    next = bottomLeftVertex;
    do
    {
        next = boundaryConnections[next];
        //Check whether the boundary coordinate is clamped to the upper edge
        if (textureCoordinates[next][1] != 0)
        {
            //Redo the barycentric mapping, ensuring it's to the correct edge
            Cartesian3 vertex = vertices[next];

            float distToLeft;
            if(axesOrientation == ZY_COORDINATES)
                distToLeft = std::abs(vertex.z - lowestZ);
            else
                distToLeft = std::abs(vertex.x - lowestX);

            float zCoordinate = distToLeft / squareWidth;

            //Update the coordinate
            textureCoordinates[next] = {zCoordinate, 0, 0};
        }
    } while(next != bottomRightVertex);


    //The boundary vertices are now sorted, so we apply the interior vertices
    for(auto pair : interiorVertices)
        textureCoordinates[pair.first] = {0.5, 0.5, 0};

    textureCoords = textureCoordinates;
}

//Returns the number of neighbours a vertex has
int AttributedObject::neighbours(int vertexID)
{
    int neighbourCount = 0;

    int firstEdge = firstDirectedEdge[vertexID];
    int currentEdge = firstEdge;
    do
    {
        //Increment the neighbour counter
        neighbourCount++;

        //Get the other half of the edge
        int other = otherHalf[currentEdge];

        //Get face and face index for other half edge
        int edgeFace = other / 3;
        int edgeFaceIndex = other % 3;

        //Go to the next edge
        currentEdge = (edgeFace * 3) + ((edgeFaceIndex + 1) % 3);

    //Keep repeating until we're back at the first edge (completed the cycle)
    } while (currentEdge != firstEdge);

    return neighbourCount;
}

void AttributedObject::performFloater()
{

    //Before performing the algorithm we need to know the boundary and interior vertices
    findBoundary();

    //We use the interior and boundary vertices to create patches
    //createPatches();
int numLoops = 0;
    //Cycle around the boundary vertices
    {
        auto test = boundaryConnections;

        int start = boundaryConnections.begin()->first;
        int current = start;



        int next = boundaryConnections[current];
        do
        {

            //std::cout << current << " -> " << next << std::endl;

            current = next;
            next = boundaryConnections[current];
            test[current] = -1;
            numLoops++;

        } while (start != next);

        int a = 5;
    }
    //Using these values, we can assign the initial texture coordinates
    mapBoundaryVertices();

    //We need two vectors, the old coordinates and the updated coordinates
    std::vector<Cartesian3> oldCoordinates, updatedCoordinates;

    //Set oldCoordinates to be the initial textureCoordinates
    oldCoordinates = textureCoords;
    updatedCoordinates = oldCoordinates;

    bool converged = true;

    //We repeat the algorithm until all interior vertices have converged
    do
    {

        //Loop through interior vertices
        for(auto pair : interiorVertices)
        {
            //Initialise a_ij
            Cartesian3 newPoint = {0,0,0};

            //Get the vertexID from interiorExterior
            int vertexID = pair.first;

            //Get the neighbours of that vertex
            int neighbourCount = neighbours(vertexID);

            //Uniform weighting
            float weight = 1.0 / neighbourCount;

            //Loop through the neighbours of the interior vertex
            //Add the position of the weighting + the weighting
            int firstEdge = firstDirectedEdge[vertexID];
            int currentEdge = firstEdge;

            do
            {
                //Get the current vertex
                int face = currentEdge / 3;
                int localIndex = currentEdge % 3;

                //Get the vertexID
                int newVertexID = faceVertices[(3*face) + localIndex];

                //Add the texture coordinate to the point
                newPoint = newPoint + (weight * oldCoordinates[newVertexID]);

                //We've added the neighbour, now we get the next neighbour
                int other = otherHalf[currentEdge];

                //Get face and local index to go to the next edge
                int otherFace = other / 3;
                int otherLocalIndex = other % 3;

                currentEdge = (otherFace * 3) + ((otherLocalIndex + 1) % 3);

            } while (firstEdge != currentEdge);

            //After adding all neighbours, we take the new point and store it in the updatedCoordinates
            updatedCoordinates[vertexID] = newPoint;
        }

        //At the conclusion, we test for convergence - see whether all interior vertices have found their places

        float tolerance = 0.001;

        for (auto pair : interiorVertices)
        {
            Cartesian3 difference = updatedCoordinates[pair.first] - oldCoordinates[pair.first];

            bool changed = false;
            //Check whether the point has moved greater than the tolerance
            if (difference.length() >= tolerance)
            {
                changed = true;
                converged = false;

                break;
            }

            //If there has been no change, then all values are below tolerance
            if (changed == false)
                converged = true;
        }

        //Next iteration, we set the old coordinates to be the updated coordinates
        oldCoordinates = updatedCoordinates;

    } while (converged == false);


    //Set textureCoords to be the most up to date texture coordinates
    textureCoords = updatedCoordinates;
}

void AttributedObject::calculateNormals()
{
    //Resize normals to fit all vertices
    normals.resize(vertices.size(), {0,0,0});

    //Loop through all faces
    for(int i = 0; i < faceVertices.size(); i+= 3)
    {
        int vertex1ID = faceVertices[i];
        int vertex2ID = faceVertices[i+1];
        int vertex3ID = faceVertices[i+2];

        Cartesian3 vertex1 = vertices[vertex1ID];
        Cartesian3 vertex2 = vertices[vertex2ID];
        Cartesian3 vertex3 = vertices[vertex3ID];

        //Calculate the face normal
        Cartesian3 edge1 = vertex2 - vertex1;
        Cartesian3 edge2 = vertex3 - vertex1;

        Cartesian3 faceNormal = edge1.cross(edge2).unit();

        //Add the face normal to each vertex
        normals.at(vertex1ID) = normals.at(vertex1ID) + faceNormal;
        normals.at(vertex2ID) = normals.at(vertex2ID) + faceNormal;
        normals.at(vertex3ID) = normals.at(vertex3ID) + faceNormal;

    }

    //After every face normal is added, we normalise every value to get vertex normals
    for(int i = 0; i < normals.size(); i++)
    {
        if(! (normals.at(i).x == 0 && normals.at(i).y == 0 && normals.at(i).z == 0))
           {
                normals.at(i) = normals.at(i).unit();
           }

        //Convert from [-1, 1] to [0,1]
        normals.at(i) = (normals.at(i) + Cartesian3(1,1,1)) * 0.5;

    }

}
