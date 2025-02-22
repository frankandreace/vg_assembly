# bdsg import
import json
from bdsg.bdsg import SnarlDistanceIndex
from bdsg.bdsg import PackedGraph
import assembler.reparser as relp
from assembler.gaf_reader import GafReader

# package import
from assembler.constants import *

from assembler.node import Node
from assembler.anchor import Anchor

# other imports
import time
from sys import stderr

class AnchorBuilder:
    """
    This class produces a Dictionary containing anchors in the pangenome graph.
    Anchors are small paths derived from bubbles in the graph (snarls) and are used by Shasta to
    phase and assembly the reads in a sample.
    A path is a succession of nodes and orientations in a bidirected graph, as described in
    https://github.com/lh3/gfatools/blob/master/doc/rGFA.md.
    The anchor dictionary has as a key a "sentinel" node in the anchor path and as value a list of tuples. In each tuple in the first position there is an anchor having the that node as sentinel and in the second position an empty list that will contain the reads associated to the anchor.
    The anchor is a list of node_handles as in the packedgraph implementation of the bdsg library (https://bdsg.readthedocs.io/en/master/index.html).
    NOTE: THIS IS GOING TO CHANGE AS I WON"T STORE ANYMORE THE HANDLES but only their useful info (node_id, orientation, node_length) and the list of reads associated to the anchors will be stored in another dictionary.

    Path example [handle_####,handle_####,handle####] -> >1>2>3 with [>,< ==  node orientation][node_id]
    Anchors example:>1>2>3      Node with ID == 2 is the sentinel of this path
                    >1>2>3>4    Node with ID == 2 is the sentinel of this path
                    <4<3<2<1    Node with ID == 2 is the sentinel of this path
    """

    def __init__(self, alignment_gaf) -> None:
        # useful initialization objects
        self.graph = PackedGraph()
        self.index = SnarlDistanceIndex()

        # important generated_data
        self.leaf_snarls: list = []
        self.leafsnarl_boundaries_dict: dict = {}   # dictionary with leaf snarl net_handle mapped to tuple with snarl boundaries (start_node_id, end_node_id)
        self.leafsnarl_boundary_nodes_inside_dict: dict = {}
        self.snarl_boundaries: list = [dict(), dict()]
        self.anchors_to_reads: dict = {}            # Dictionary with anchors and list of all reads containing the anchor (not necessarily bp matched)
        self.anchors_to_bpmatched_reads: dict = {}  # Dictionary with anchors and list of all bp-matched reads

        # temporary variables to store data between functions
        self.contains_child_snarls: bool = False
        # self.keep_path_scan: bool = True
        # self.path_orientation: bool = True
        self.current_snarl_start: int = 0
        # self.peek_orientations = []
        # self.count_in_path: bool = True
        self.snarl_id: int = 0

        self.current_anchor: Anchor = Anchor()
        # self.curr_path_name = ""
        # self.verbose = False
        self.ref_path_name = "CHM13".casefold()
        # self.path_names = []

        # variables used for debugging
        self.num_usable_bubbles = 0
        self.num_used_bubbles = 0

        # alignment related
        self.gaf_reader = GafReader(alignment_gaf)

    def build(self, packed_graph_path: str, index_path: str) -> None:
        """
        Deserializes the packedGraph and SnarlIndexes generates using vg. Does not return anything

        Parameters
        ----------
        packed_graph_path : string
            Path to packedGraph object (.vg)
        index_path : string
            Path to SnarlIndex object (.dist)
        """
        t0=time.time()
        self.graph.deserialize(packed_graph_path)
        self.index.deserialize(index_path)
        print(f"Graph files loaded in {time.time()-t0:.2f}", file=stderr)

    
    def get_leaf_snarls_with_boundaries(self) -> None:
        """
        It identifies all leaf snarls in the graph/subgraph by traversing snarlDistanceIndex decomposition. It also stores snarl boundary nodes.

        Returns
        -------
        None
        """

        # Compute leaf snarls
        t0 = time.time()
        if len(self.leaf_snarls) == 0:
            self.process_snarls()
        print(
            f"Leaf Snarls Computed in {time.time()-t0:.2f}s",
            file=stderr,
        )

        # compute boundaries of leaf snarls
        t1 = time.time()
        self.create_snarl_boundaries_dict()
        print(
            f"Snarl Boundaries computed in {time.time()-t1:.2f}s",
            file=stderr,
        )


    def process_snarls(self) -> None:
        """
        This function traverses the whole Snarl Tree index and stores the leaf snarls into a list for future processing into anchors.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.index.traverse_decomposition(
            self.check_leaf_snarl_iteratee,  # snarl_iteratee
            lambda x: True,  #  chain_iteratee
            lambda y: True,  # node_iteratee
        )
        return None

    
    def check_leaf_snarl_iteratee(self, net_handle) -> bool:
        """
        This function is called on the snarl tree traversal (process_snarls function) when the pointer is on a snarl net_handle. It verifies if the snarl is a leaf snalr (does not contain inside it another snarl like a matrioska).
        If the snarl is found to be a leaf snarl, it appends its net_handle to a list of valid snarl handles to then process them to generate anchors.
        It returns True to keep the iteration going and do not stop it.

        Parameters
        ----------
        net_handle : object
            net_handle object from a SnarlIndex

        Returns
        -------
        bool
        True as iteration has to continue
        """

        self.contains_child_snarls = False

        snarl_children: list = []
        self.index.for_each_child(
            net_handle, lambda y: snarl_children.append(y) or True
        )

        for s_c in snarl_children:
            self.index.for_each_child(s_c, self.check_snarl_in_children_iteratee)

        if not self.contains_child_snarls:
            self.leaf_snarls.append(net_handle)

        return True


    def check_snarl_in_children_iteratee(self, child_net_handle) -> bool:
        """
        It iterates on the children of each snarl child (check ) to verify that the snalr does not contain any other snarl and is therefore a snarl leave. If it seesa a snarl it sets the variable contains_child_snarls as true.
        The return True/ False parameter is used to continue or stop the iteration calling this function. It stops if it finds a snarl else continue to check the handles of the snarl childern.

        Parameters
        ----------
        child_net_handle : object
            net_handle object from a SnarlIndex

        Returns
        -------
        bool
        True if iteration has to continue else False
        """

        if self.index.is_snarl(child_net_handle):
            self.contains_child_snarls = True
            return False

        return True


    def create_snarl_boundaries_dict(self):

        for _, snarl_net_handle in enumerate(self.leaf_snarls):
            snarl_boundary, nodes_in_snarl = self.get_snarl_boundaries_handle_chuncked_graph(snarl_net_handle)
            self.leafsnarl_boundaries_dict[snarl_net_handle] = snarl_boundary
            # Save a dictionary with {snarl_boundary: nodes in snarl}
            self.leafsnarl_boundary_nodes_inside_dict[snarl_boundary] = nodes_in_snarl

    
    # Returns a tuple containing the start and end boundaries of a snarl
    def get_snarl_boundaries_handle_chuncked_graph(self, snarl_net_handle) -> tuple:
        """
        This function takes a snarl net_handle and returns the boundary nodes of the snarl, i.e. preceding and succeeding the snarl. 
        It also creates the self.snarl_boundaries list of 2 dictionaries (for FORWARD and REVERSE)

        Parameters
        ----------
        snarl_net_handle: obj
        net_handle object of the snarl

        Returns
        -------
        snarl_boundary : tuple
        the node_id of the nodes preceding and succeding the snarl
        nodes_in_snarl: list
        list of nodes (node_ids) inside the snarl
        """

        # get start and end boundary of a snarl
        start_bound = self.index.get_start_bound(snarl_net_handle)
        end_bound = self.index.get_end_bound(snarl_net_handle)
        # Convert boundary net handles to node handles
        start_node_handle = self.index.get_handle(start_bound, self.graph)
        end_node_handle = self.index.get_handle(end_bound, self.graph)

        # Now create a FORWARD and REVERSE dictionary
        start_id = self.graph.get_id(start_node_handle)
        end_id = self.graph.get_id(end_node_handle)
        snarl_boundary = (min(start_id, end_id), max(start_id, end_id))

        self.snarl_boundaries[FORWARD_DICTIONARY][snarl_boundary[0]] = (
            snarl_boundary[1],
            snarl_net_handle
        )
        self.snarl_boundaries[REVERSE_DICTIONARY][snarl_boundary[1]] = (
            snarl_boundary[0],
            snarl_net_handle
        )

        # Nodes within the snarl boundary (doesn't include boundary nodes)
        nodes_in_snarl = self.get_nodes_in_snarl(snarl_net_handle)
        
        return snarl_boundary, nodes_in_snarl


    def get_nodes_in_snarl(self, snarl_net_handle) -> list:
        nodes_inside = []

        self.index.traverse_decomposition_helper(
            snarl_net_handle,  
            snarl_iteratee=lambda s: True,  # Ignore snarls
            chain_iteratee=lambda c: True,  # Ignore chains
            node_iteratee=lambda n: nodes_inside.append(self.graph.get_id(self.index.get_handle(n, self.graph))) or True
        )

        return nodes_inside  # Return the collected node IDs


    def dump_forward_reverse_snarl_boundaries_dict(self, output_prefix):
        print(f"Printing to {output_prefix}.forward_dict.csv")
        with open(f"{output_prefix}.forward_dict.csv", "w") as f:
            f.write("snarl_start,snarl_end,snarl_net_handle\n")
            for el in self.snarl_boundaries[FORWARD_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[FORWARD_DICTIONARY][el][0]},{self.snarl_boundaries[FORWARD_DICTIONARY][el][1]}", file=f
                )
        print(f"Printing to {output_prefix}.reverse_dict.csv")
        with open(f"{output_prefix}.reverse_dict.csv", "w") as f:
            f.write("snarl_start,snarl_end,snarl_net_handle\n")
            for el in self.snarl_boundaries[REVERSE_DICTIONARY]:
                print(
                    f"{el},{self.snarl_boundaries[REVERSE_DICTIONARY][el][0]},{self.snarl_boundaries[REVERSE_DICTIONARY][el][1]}", file=f
                )


    def dump_snarl_dictionary(self, file_path):
        print(f"Dumping snarl boundary to nodes inside snarl dictionary to TSV file")
        with open(file_path, "w") as f:
            f.write("snarl_start\tsnarl_end\tnodes_in_snarl\n")
            for snarl_boundary, nodes_in_snarl in self.leafsnarl_boundary_nodes_inside_dict.items():
                print(
                    f"{snarl_boundary[0]}\t{snarl_boundary[1]}\t{nodes_in_snarl}", file=f
                )
    

    ############################################################################################################################################################
    ## Alignment processing....
    def process(self, debug_outfile, reads_out_file):
        """
        It reads the gaf file line by line and if the line is valid, it processes it to find anchors that align to it.
        """
        times = []
        with open(debug_outfile, 'w') as debug:
            print("READ_ID\tANCHOR\tIS_MATCHING_NODES\tIS_BASELEVEL_ALIGNED", file=debug)
            with open(reads_out_file, 'w') as reads_f:
                print("READ_ID\tREAD_LEN\tRELATIVE_STRAND\tPATH_START\tPATH_END\tNODES_LIST\tORIENTATION_LIST\tCS_LINE", file=reads_f)
            for line in self.gaf_reader.get_lines():
                t0 = time.time()
                parsed_data = relp.processGafLine(line, reads_out_file)
                if parsed_data:
                    print(
                        f"PROCESSING READ {parsed_data[0]} ...",
                        end="\n",
                        flush=True,
                        file=stderr,
                    )
                    self.processGafLine(parsed_data, debug)
                    t1 = time.time()
                    print(f"Done in {t1-t0}.", file=stderr)
                    times.append(t1-t0)


    def processGafLine(self, alignment_l: list, debug_file):

        # walk through the nodes keeping the position in the node list
        # and the total length of the nodes
        read_id = alignment_l[READ_POSITION]
        forward_dict = self.snarl_boundaries[FORWARD_DICTIONARY]
        reverse_dict = self.snarl_boundaries[REVERSE_DICTIONARY]

        # # Verifying that the nodes coming from the alingment are in the graph I am using
        # if not self.graph.has_node(node_id):
        #     continue
        
        # Traversing each node in read alignment path
        for position, node_id in enumerate(alignment_l[NODE_POSITION]):
            walked_in_read_bps = 0
            walked_inside_anchor_bps = 0
            got_anchor = False

            # Check relative orientation of the read
            if alignment_l[STRAND_POSITION]:   # If forward strand
                if node_id in forward_dict:
                    anchor_start_node = node_id
                    anchor_end_node = forward_dict[anchor_start_node][0]
                    anchor_start_pos_in_read = anchor_end_pos_in_read = position
                    # Nodes inside the snarl
                    nodes_inside = self.leafsnarl_boundary_nodes_inside_dict[(anchor_start_node,anchor_end_node)]
                    for i, next_node in enumerate(alignment_l[NODE_POSITION][position:position+len(nodes_inside)+1]):
                        if next_node == anchor_end_node:
                            anchor_end_pos_in_read = anchor_start_pos_in_read + i
                            node_handle = self.graph.get_handle(node_id)
                            length = self.graph.get_length(node_handle)
                            walked_inside_anchor_bps+=length
                            ## Here also check for base level alignments
                            cigar_string = self.get_cigar(walked_in_read_bps, alignment_l[CIGAR_POSITION], read_id, walked_inside_anchor_bps)
                            got_anchor=True
                            break
                        node_handle = self.graph.get_handle(node_id)
                        length = self.graph.get_length(node_handle)
                        walked_inside_anchor_bps+=length

                    if got_anchor:
                        # Save anchor with read info in a dictionary
                        anchor = ">".join(map(str, alignment_l[NODE_POSITION][anchor_start_pos_in_read:anchor_end_pos_in_read+1]))
                        if anchor not in self.anchors_to_reads:
                            self.anchors_to_reads[anchor] = []
                        self.anchors_to_reads[anchor].append([read_id, alignment_l[STRAND_POSITION], anchor_start_pos_in_read, anchor_end_pos_in_read, walked_inside_anchor_bps, cigar_string])
            
                node_handle = self.graph.get_handle(node_id)
                length = self.graph.get_length(node_handle)
                walked_in_read_bps+=length
            
            else:
                if node_id in reverse_dict:
                    anchor_start_node = node_id
                    anchor_end_node = reverse_dict[anchor_start_node][0]
                    anchor_start_pos_in_read = anchor_end_pos_in_read = position
                    # Nodes inside the snarl
                    nodes_inside = self.leafsnarl_boundary_nodes_inside_dict[(anchor_end_node,anchor_start_node)]
                    for i, next_node in enumerate(alignment_l[NODE_POSITION][position:position+len(nodes_inside)+1]):
                        if next_node == anchor_end_node:
                            anchor_end_pos_in_read = anchor_start_pos_in_read + i
                            node_handle = self.graph.get_handle(node_id)
                            length = self.graph.get_length(node_handle)
                            walked_inside_anchor_bps+=length
                            ## Here also check for base level alignments
                            cigar_string = self.get_cigar(walked_in_read_bps, alignment_l[CIGAR_POSITION], read_id, walked_inside_anchor_bps)
                            got_anchor=True
                            break
                        node_handle = self.graph.get_handle(node_id)
                        length = self.graph.get_length(node_handle)
                        walked_inside_anchor_bps+=length

                    if got_anchor:
                        # Save anchor with read info in a dictionary
                        anchor = ">".join(map(str, alignment_l[NODE_POSITION][anchor_start_pos_in_read:anchor_end_pos_in_read+1][::-1]))
                        if anchor not in self.anchors_to_reads:
                            self.anchors_to_reads[anchor] = []
                        self.anchors_to_reads[anchor].append([read_id, alignment_l[STRAND_POSITION], anchor_start_pos_in_read, anchor_end_pos_in_read, walked_inside_anchor_bps, cigar_string])
                
                node_handle = self.graph.get_handle(node_id)
                length = self.graph.get_length(node_handle)
                walked_in_read_bps+=length



    def get_cigar(self, walked_in_read_bps: int, cs_line: list, read_id: str, walked_inside_anchor_bps: int):
        """
        It uses the parsed cs tag from the gaf to verify that the anchor and the path match at the sequence level.
        """

        # Reconstruct the read sequence from the CS tag
        # We know we have walked `walked_in_read_bps` bases in the read, so use that to extract 

        anchor_end_bps = walked_in_read_bps + walked_inside_anchor_bps
        cigar_string = ""
        for symbol, bps in cs_line:
            walked_in_read_bps -= bps
            anchor_end_bps -= bps
            if walked_in_read_bps < 0:
                start_index_inside_cs_tag = max(0, walked_in_read_bps+bps)
                end_index_inside_cs_tag = min(bps, anchor_end_bps+bps)
                cigar_string += (end_index_inside_cs_tag-start_index_inside_cs_tag)*symbol
                if anchor_end_bps < 0:
                    break
        
        return(cigar_string)
    

    def dump_bp_matched_reads(self):
        """
        Find reads with exact base level matches to the anchors in the graph.
        """
        
        for anchor, reads in self.anchors_to_reads.items():
            for read in reads:
                cigar = read[5]
                if set(cigar) == {":"}:   # perfect base-level matches
                    if anchor not in self.anchors_to_bpmatched_reads:
                        self.anchors_to_bpmatched_reads[anchor] = []
                    self.anchors_to_bpmatched_reads[anchor].append(read[0:4])


    def dump_anchors_to_json_for_shasta(self, anchors_json_file: str):
        """
        Dumps anchors to a JSON file as required by Shasta
        """
        
        json_list = []
        for anchor, reads in self.anchors_to_bpmatched_reads.items():
            reads.insert(0, anchor)
            json_list.append(reads)

        with open(anchors_json_file, "w", encoding="utf-8") as f:
            json.dump(json_list, f, ensure_ascii=False)


    ### PRINTING FUNCTIONS FOR DEBUG - VISUALIZATION ###

    # def print_anchors_from_dict(self, file_name) -> None:
    #     with open(file_name, "w") as f:
    #         for sentinel, anchor_list in self.sentinel_to_anchor.items():
    #             for anchor in anchor_list:
    #                 print(
    #                     f"Sentinel: {sentinel} ; Anchor : {anchor!r} ; Bandage : {anchor.bandage_representation()}",
    #                     file=f,
    #                 )

    # def print_traversal(self, traversal: list) -> None:
    #     for node in traversal:
    #         # node_handle = self.index.get_handle(node_net_handle, self.graph)
    #         direction = ">" if node.orientation == True else "<"
    #         print(f"{direction}{node.id}", end="")
    #     print()

    # def print_tree_snarl_iteratee(self, snarl_handle) -> bool:
    #     print("Snarl:", self.index.net_handle_as_string(snarl_handle), flush=True)
    #     return True

    # def print_tree_chain_iteratee(self, chain_handle) -> bool:
    #     print("Chain:", self.index.net_handle_as_string(chain_handle), flush=True)
    #     return True

    # def print_tree_node_iteratee(self, node_handle) -> bool:
    #     print("Node:", self.index.net_handle_as_string(node_handle), flush=True)
    #     return True

    # def print_tree_structure(self) -> None:
    #     print("Printing tree structure now", flush=True)
    #     self.index.traverse_decomposition(
    #         self.print_tree_snarl_iteratee,  # snarl_iteratee
    #         self.print_tree_chain_iteratee,  #  chain_iteratee
    #         self.print_tree_node_iteratee,  # node_iteratee
    #     )

    # def print_sentinels_for_bandage(self, file) -> None:
    #     sentinel_nodes_set = set()
    #     for _, anchor_list in self.sentinel_to_anchor.items():
    #         for anchor in anchor_list:
    #             for node_h in anchor:
    #                 sentinel_nodes_set.add(node_h.id)

    #     with open(file, "w") as out_f:
    #         print("Node,color", file=out_f)
    #         for node in sentinel_nodes_set:
    #             print(f"{node},#FF0000", file=out_f)

    # def get_dict(self) -> dict:
    #     return self.sentinel_to_anchor

    # def get_path_names(self):
    #     path_names = []

    #     def collect_name(path_handle):
    #         path_names.append(self.graph.get_path_name(path_handle))
    #         return True

    #     self.graph.for_each_path_handle(collect_name)
    #     return path_names

    # def print_dict_sizes(self, out_f) -> None:
    #     with open(out_f, "w") as f:
    #         print(f"Sentinel_node\tsnarl_id\tAnchor_length\tAnchor_pos_in_ref_path\tAnchor_path\tAnchor_nodes_copypaste_bandage\tPaths_associated_with_anchor",file=f)
    #         for sentinel, anchor_list in self.sentinel_to_anchor.items():
    #             for anchor in anchor_list:
    #                 print(
    #                     f"{sentinel}\t{anchor.snarl_id}\t{anchor.baseparilength}\t{anchor.genomic_position}\t{anchor!r}\t{anchor.bandage_representation()}\t{anchor.get_reference_paths()}",
    #                     file=f,
    #                 )