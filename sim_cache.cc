#include <iostream>
#include<string>
#include<fstream>
#include <math.h>
#include<stdio.h>
using namespace std;

typedef struct{
    unsigned int t_tag, t_LRU_counter, t_dirty_bit, t_valid_bit; //Define structure for L1 cache
}table;

typedef struct{
    unsigned int* ds_tag;
}decouple_sector_tag; //Define structure for decoupled sector cache tag --- L2

typedef struct{
    unsigned int *ds_dirty_bit, *ds_valid_bit, *ds_sel_bit;
}decouple_sector_data; //Define structure for decoupled sector cache dataset --- L2

class L2CACHE
{
	public:
		unsigned int BLOCKSIZE, L2_SIZE, L2_ASSOC, sets, bo, indexsizeL2, tagsizeL2, setsL2, L2_numtag, L2_numblocks, L2En, DSEn;
	unsigned int readsL2, writesL2, read_missesL2, write_missesL2, write_backsL2, sector_miss, block_miss;
	unsigned int L2_DATA_BLOCKS, L2_ADDR_TAGS;
	table ** ArrL2; //Two dimentional array to hold the index values which will point to the blocks in that index
	decouple_sector_tag * L2_tagstore; //to hold tagstore of L2
	decouple_sector_data * L2_datastore; //to hold datastore of L2


	L2CACHE(unsigned int c_BLOCKSIZE, unsigned int c_L2_SIZE, unsigned int c_L2_ASSOC, unsigned int c_L2_DATA_BLOCKS, unsigned int c_L2_ADDR_TAGS)
	{

		BLOCKSIZE = c_BLOCKSIZE;
		L2_SIZE = c_L2_SIZE;
		L2_ASSOC = c_L2_ASSOC;
		L2_DATA_BLOCKS = c_L2_DATA_BLOCKS;
		L2_ADDR_TAGS = c_L2_ADDR_TAGS;
		if (c_L2_SIZE > 0)
		{
			L2_numtag = log2(c_L2_ADDR_TAGS);
			L2_numblocks = log2(c_L2_DATA_BLOCKS);
			sets = c_L2_SIZE / (c_BLOCKSIZE * c_L2_ASSOC * c_L2_DATA_BLOCKS); //number of sets
			bo = log2(c_BLOCKSIZE); //block offset
			indexsizeL2 = log2(sets); //number of bits in index value
			tagsizeL2 = 32 - bo - indexsizeL2 - L2_numblocks - L2_numtag;
			readsL2 = 0;
			writesL2 = 0;
			read_missesL2 = 0;
			write_missesL2 = 0;
			write_backsL2 = 0;
			sector_miss = 0;
			block_miss = 0;
			DSEn = 0;
			L2En = 0;

			if (c_L2_ADDR_TAGS > 1 || c_L2_DATA_BLOCKS > 1)
			{
				DSEn = 1;
			}
			if (c_L2_ADDR_TAGS == 1 && c_L2_DATA_BLOCKS == 1)
			{
				L2En = 1;
			}

			ArrL2 = new table * [sets]; //Define the two dimentional array
			for (unsigned int i = 0; i < sets; i++)
				ArrL2[i] = new table[L2_ASSOC];

			for (unsigned int i = 0; i < sets; i++)
			{
				for (unsigned int j = 0; j < L2_ASSOC; j++)
				{
					ArrL2[i][j].t_LRU_counter = j; //initialise the lru counters, dirty bits,valid bits and tags
					ArrL2[i][j].t_dirty_bit = 0;
					ArrL2[i][j].t_valid_bit = 0;
					ArrL2[i][j].t_tag = 0;


				}
			}

			L2_tagstore = new decouple_sector_tag[sets]; //Define the two dimentional array for L2 tag storage

			for (unsigned int i = 0; i < sets; i++)
			{
				L2_tagstore[i].ds_tag = new unsigned int[L2_ADDR_TAGS];
				for (unsigned int k = 0; k < L2_ADDR_TAGS; k++)
				{
					L2_tagstore[i].ds_tag[k] = 0;
				}
			}

			L2_datastore = new decouple_sector_data[sets]; //Define the two dimentional array for L2 data storage

			for (unsigned int i = 0; i < sets; i++)
			{
				L2_datastore[i].ds_valid_bit = new unsigned int[L2_DATA_BLOCKS];
				L2_datastore[i].ds_sel_bit = new unsigned int[L2_DATA_BLOCKS];
				L2_datastore[i].ds_dirty_bit = new unsigned int[L2_DATA_BLOCKS];
				for (unsigned int k = 0; k < L2_DATA_BLOCKS; k++)
				{
					L2_datastore[i].ds_valid_bit[k] = 0;
					L2_datastore[i].ds_dirty_bit[k] = 0;
					L2_datastore[i].ds_sel_bit[k] = 0;

				}

			}
		}

	}


	void readFromAddressL2(unsigned int);
	void writeToAddressL2(unsigned int);
	void readFromAddressDS(unsigned int);
	void writeToAddressDS(unsigned int);
	void printStatsL2();
};

class L1CACHE
{
	public:
		unsigned int BLOCKSIZE, L1_SIZE, L1_ASSOC, sets, block, bo, indexsize, tagsize;
	unsigned int reads, writes, read_misses, write_misses, write_backs, L2En, DSEn;
	L2CACHE * next; //Next level cache
	table ** Arr; //Two dimentional array to hold the index values which will point to the blocks in that index
	int l2count;

	L1CACHE(unsigned int c_BLOCKSIZE, unsigned int c_L1_SIZE, unsigned int c_L1_ASSOC, L2CACHE * c_next, unsigned int c_DSEn, unsigned int c_L2En)
	{
		BLOCKSIZE = c_BLOCKSIZE;
		L1_SIZE = c_L1_SIZE;
		L1_ASSOC = c_L1_ASSOC;
		block = c_L1_SIZE / c_BLOCKSIZE; //number of blocks
		sets = c_L1_SIZE / (c_BLOCKSIZE * c_L1_ASSOC); //number of sets
		bo = log2(c_BLOCKSIZE); //block offset
		indexsize = log2(sets); //number of bits in index value
		tagsize = 32 - bo - indexsize; //tag length
		next = c_next; //Next level(NULL for Project:1A)                                
		reads = 0;
		writes = 0;
		read_misses = 0;
		write_misses = 0;
		write_backs = 0;
		L2En = c_L2En;
		DSEn = c_DSEn;
		Arr = new table * [sets]; //Define the two dimentional array
		for (unsigned int i = 0; i < sets; i++)
			Arr[i] = new table[L1_ASSOC];

		for (unsigned int i = 0; i < sets; i++)
		{
			for (unsigned int j = 0; j < L1_ASSOC; j++)
			{
				Arr[i][j].t_LRU_counter = j; //initialise the lru counters, dirty bits,valid bits and tags
				Arr[i][j].t_dirty_bit = 0;
				Arr[i][j].t_valid_bit = 0;
				Arr[i][j].t_tag = 0;
			}
		}
	}
	void readFromAddressL1(unsigned int);
	void writeToAddressL1(unsigned int);
	void printStatsL1();
};


void L2CACHE::readFromAddressDS(unsigned int addr)
{
	readsL2++;
	unsigned int tag = addr >> (bo + indexsizeL2 + L2_numtag + L2_numblocks); //extract the tag from input address
	unsigned int index = ((1 << indexsizeL2) - 1) & (addr >> (bo + L2_numblocks)); //extract index
	unsigned int Ntag = ((1 << L2_numtag) - 1) & (addr >> (bo + L2_numblocks + indexsizeL2)); //extract N
	unsigned int Pblock = ((1 << L2_numblocks) - 1) & (addr >> bo);

	int r_hit = 0;
	if (L2_tagstore[index].ds_tag[Ntag] == tag && L2_datastore[index].ds_valid_bit[Pblock] == 1 && L2_datastore[index].ds_sel_bit[Pblock] == Ntag) //find if hit or miss and if hit find the hit block
	{
		r_hit = 1;
	}

	if (r_hit == 0) //miss case
	{

		if (L2_datastore[index].ds_dirty_bit[Pblock] == 1)
		{
			write_backsL2++;
			L2_datastore[index].ds_dirty_bit[Pblock] = 0;
		}

		unsigned int FoundValid = 0;
		for (unsigned int i = 0; i < L2_DATA_BLOCKS; i++)
		{
			if (L2_datastore[index].ds_valid_bit[i] == 1)
				FoundValid = 1;
		}

		if (FoundValid == 1)
			block_miss++;
		else
			sector_miss++;

		if (L2_tagstore[index].ds_tag[Ntag] != tag)
		{
			for (unsigned int i = 0; i < L2_DATA_BLOCKS; i++)
			{
				if (L2_datastore[index].ds_sel_bit[i] == Ntag)
				{
					if (L2_datastore[index].ds_dirty_bit[i] == 1)
					{
						write_backsL2++;
						L2_datastore[index].ds_dirty_bit[i] = 0;
					}
					L2_datastore[index].ds_valid_bit[i] = 0;
					L2_datastore[index].ds_sel_bit[i] = 0;
				}

			}
			L2_tagstore[index].ds_tag[Ntag] = tag;
			L2_datastore[index].ds_valid_bit[Pblock] = 1;
			L2_datastore[index].ds_sel_bit[Pblock] = Ntag;
			L2_datastore[index].ds_dirty_bit[Pblock] = 0;

		}

		else if (L2_tagstore[index].ds_tag[Ntag] == tag)
		{
			L2_datastore[index].ds_valid_bit[Pblock] = 1;
			L2_datastore[index].ds_sel_bit[Pblock] = Ntag;
			L2_datastore[index].ds_dirty_bit[Pblock] = 0;
		}
		read_missesL2++;

	}

}

void L2CACHE::writeToAddressDS(unsigned int z)
{
	writesL2++;
	unsigned int tag = z >> (bo + indexsizeL2 + L2_numtag + L2_numblocks); //extract the tag from input address
	unsigned int index = ((1 << indexsizeL2) - 1) & (z >> (bo + L2_numblocks)); //extract index
	unsigned int Ntag = ((1 << L2_numtag) - 1) & (z >> (bo + L2_numblocks + indexsizeL2)); //extract N
	unsigned int Pblock = ((1 << L2_numblocks) - 1) & (z >> bo);
	int r_hit = 0;
	if (L2_tagstore[index].ds_tag[Ntag] == tag && L2_datastore[index].ds_valid_bit[Pblock] == 1 && L2_datastore[index].ds_sel_bit[Pblock] == Ntag) //find if hit or miss and if hit find the hit block
	{
		r_hit = 1;
		L2_datastore[index].ds_dirty_bit[Pblock] = 1;
	}

	if (r_hit == 0) //miss case
	{


		if (L2_datastore[index].ds_dirty_bit[Pblock] == 1)
		{
			write_backsL2++;
			L2_datastore[index].ds_dirty_bit[Pblock] = 0;
		}

		unsigned int FoundValid = 0;
		for (unsigned int i = 0; i < L2_DATA_BLOCKS; i++)
		{
			if (L2_datastore[index].ds_valid_bit[i] == 1)
				FoundValid = 1;
		}

		if (FoundValid == 1)
			block_miss++;
		else
			sector_miss++;


		if (L2_tagstore[index].ds_tag[Ntag] != tag)
		{
			for (unsigned int i = 0; i < L2_DATA_BLOCKS; i++)
			{
				if (L2_datastore[index].ds_sel_bit[i] == Ntag)
				{
					if (L2_datastore[index].ds_dirty_bit[i] == 1)
					{
						write_backsL2++;
						L2_datastore[index].ds_dirty_bit[i] = 0;
					}
					L2_datastore[index].ds_valid_bit[i] = 0;
					L2_datastore[index].ds_sel_bit[i] = 0;
				}

			}
			L2_tagstore[index].ds_tag[Ntag] = tag;
			L2_datastore[index].ds_valid_bit[Pblock] = 1;
			L2_datastore[index].ds_sel_bit[Pblock] = Ntag;
			L2_datastore[index].ds_dirty_bit[Pblock] = 1;

		}

		else if (L2_tagstore[index].ds_tag[Ntag] == tag)
		{
			L2_datastore[index].ds_valid_bit[Pblock] = 1;
			L2_datastore[index].ds_sel_bit[Pblock] = Ntag;
			L2_datastore[index].ds_dirty_bit[Pblock] = 1;
		}
		write_missesL2++;
	}
}


void L2CACHE::readFromAddressL2(unsigned int addr)
{
	readsL2++;
	int r_hit_block = 0;
	unsigned int tag = addr >> (bo + indexsizeL2); //extract the tag from input address
	unsigned int index = ((1 << indexsizeL2) - 1) & (addr >> (bo));

	int r_hit = 0;
	for (unsigned int j = 0; j < L2_ASSOC; j++)
	{
		if (ArrL2[index][j].t_tag == tag) //find if hit or miss and if hit find the hit block
		{
			r_hit = 1;
			r_hit_block = j;
		}
	}

	if (r_hit == 1) //hit case
	{
		unsigned int oldLRU;
		oldLRU = ArrL2[index][r_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L2_ASSOC; i++)
		{
			if (ArrL2[index][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
				ArrL2[index][i].t_LRU_counter += 1; //is less than hit block counter value
		}
		ArrL2[index][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
	}


	if (r_hit == 0) //miss case
	{
		read_missesL2 = read_missesL2 + 1;

		int invalid = 0;
		for (unsigned int vl = 0; vl < L2_ASSOC; vl++)
		{
			if (ArrL2[index][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;

			}
		}


		if (invalid == 1) //case when only one block is empty
		{
			for (unsigned int i = 0; i < L2_ASSOC; i++)
			{
				if (ArrL2[index][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit 
				{
					for (unsigned int j = 0; j < L2_ASSOC; j++)
					{

						ArrL2[index][j].t_LRU_counter++;
					}

					ArrL2[index][i].t_tag = tag;
					ArrL2[index][i].t_LRU_counter = 0;
					ArrL2[index][i].t_valid_bit = 1;
				}
			}

		}

		if (invalid == 0) //when all are filled, check highest LRU_count and update tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = ArrL2[index][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L2_ASSOC; j++)
			{
				if (ArrL2[index][j].t_LRU_counter > temp)
				{
					temp = ArrL2[index][j].t_LRU_counter;
					highest_LRU = j;
				}
			}

			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{

				ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag;
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			if (ArrL2[index][highest_LRU].t_dirty_bit == 1)
				write_backsL2++;
			ArrL2[index][highest_LRU].t_dirty_bit = 0;

		}

		if (invalid != 0 && invalid != 1) //case when multiple blocks are empty
		{
			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_valid_bit == 0 && ArrL2[index][k].t_LRU_counter > currentHighest) //check which are empty and find the
				{ //one with highest lru count
					currentHighest = ArrL2[index][k].t_LRU_counter;
					highest_LRU = k;
				}

			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_LRU_counter < ArrL2[index][highest_LRU].t_LRU_counter)
					ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag; //update tag,lru counter and valid bit 
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			ArrL2[index][highest_LRU].t_valid_bit = 1;
		}

	}
}





void L2CACHE::writeToAddressL2(unsigned int z)
{
	writesL2++;
	int w_hit_block = 0;
	unsigned int tag = z >> (bo + indexsizeL2); //extract tag from address
	unsigned int index = ((1 << indexsizeL2) - 1) & (z >> (bo)); //extract index from address
	//unsigned int block_bits =  ((1<<bo)-1) & (z>>1); 
	int w_hit = 0;
	for (unsigned int j = 0; j < L2_ASSOC; j++)
	{
		if (ArrL2[index][j].t_tag == tag && ArrL2[index][j].t_valid_bit == 1) //find if hit or miss and if hit find the hit block
		{
			w_hit = 1;
			w_hit_block = j;
		}
	}

	if (w_hit == 1) //hit case
	{
		unsigned int oldLRU;
		oldLRU = ArrL2[index][w_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L2_ASSOC; i++)
		{
			if (ArrL2[index][i].t_LRU_counter < oldLRU) //update counters of other blocks which have lru counter value
				ArrL2[index][i].t_LRU_counter += 1; //less than hit block counter value
		}
		ArrL2[index][w_hit_block].t_LRU_counter = 0; //put lru counter value of hit block to zero and put dirty flag to 1
		ArrL2[index][w_hit_block].t_dirty_bit = 1;
	}


	if (w_hit == 0) //miss case
	{
		write_missesL2++;

		unsigned int invalid = 0;
		for (unsigned int vl = 0; vl < L2_ASSOC; vl++)
		{
			if (ArrL2[index][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;

			}
		}


		if (invalid == 1) //case when only one block is empty
		{
			for (unsigned int i = 0; i < L2_ASSOC; i++)
			{
				if (ArrL2[index][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit when only one block is invalid(empty)
				{
					for (unsigned int j = 0; j < L2_ASSOC; j++)
					{

						ArrL2[index][j].t_LRU_counter++;
					}

					ArrL2[index][i].t_tag = tag;
					ArrL2[index][i].t_LRU_counter = 0;
					ArrL2[index][i].t_valid_bit = 1;
					ArrL2[index][i].t_dirty_bit = 1;


				}
			}
		}

		if (invalid == 0) //when none are empty, check highest LRU_count and add tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = ArrL2[index][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L2_ASSOC; j++)
			{
				if (ArrL2[index][j].t_LRU_counter > temp)
				{
					temp = ArrL2[index][j].t_LRU_counter;
					highest_LRU = j;
				}
			}

			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{

				ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag;
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			if (ArrL2[index][highest_LRU].t_dirty_bit == 1)
				write_backsL2++;
			ArrL2[index][highest_LRU].t_dirty_bit = 1;


		}

		if (invalid != 0 && invalid != 1) //when multiple are empty 
		{
			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_valid_bit == 0 && ArrL2[index][k].t_LRU_counter > currentHighest) //find the empty blocks and find the
				{ //one with highest lru counter value
					currentHighest = ArrL2[index][k].t_LRU_counter;
					highest_LRU = k;
				}

			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L2_ASSOC; k++)
			{
				if (ArrL2[index][k].t_LRU_counter < ArrL2[index][highest_LRU].t_LRU_counter)
					ArrL2[index][k].t_LRU_counter++;
			}

			ArrL2[index][highest_LRU].t_tag = tag;
			ArrL2[index][highest_LRU].t_LRU_counter = 0;
			ArrL2[index][highest_LRU].t_valid_bit = 1;
			ArrL2[index][highest_LRU].t_dirty_bit = 1;
		}

	}

}


void L1CACHE::readFromAddressL1(unsigned int addr)
{
	reads++;
	int r_hit_block = 0;
	unsigned int tagr = addr >> (bo + indexsize); //extract the tag from input address
	unsigned int indexr = ((1 << indexsize) - 1) & (addr >> (bo));

	int r_hit = 0;
	for (unsigned int j = 0; j < L1_ASSOC; j++)
	{
		if (Arr[indexr][j].t_tag == tagr) //find if hit or miss and if hit find the hit block
		{
			r_hit = 1;
			r_hit_block = j;
		}
	}

	if (r_hit == 1) //hit case
	{
		unsigned int oldLRU;
		oldLRU = Arr[indexr][r_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L1_ASSOC; i++)
		{
			if (Arr[indexr][i].t_LRU_counter < oldLRU) //update lru counters if the value of counters
				Arr[indexr][i].t_LRU_counter += 1; //is less than hit block counter value
		}
		Arr[indexr][r_hit_block].t_LRU_counter = 0; //put lru counter of hit block to zero
	}


	if (r_hit == 0) //miss case
	{
		read_misses++;
		// next->readFromAddressL2(addr);
		int invalid = 0;
		for (unsigned int vl = 0; vl < L1_ASSOC; vl++)
		{
			if (Arr[indexr][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;

			}
		}


		if (invalid == 1) //case when only one block is empty
		{
			for (unsigned int i = 0; i < L1_ASSOC; i++)
			{
				if (Arr[indexr][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit 
				{
					for (unsigned int j = 0; j < L1_ASSOC; j++)
					{

						Arr[indexr][j].t_LRU_counter++;
					}

					Arr[indexr][i].t_tag = tagr;
					Arr[indexr][i].t_LRU_counter = 0;
					Arr[indexr][i].t_valid_bit = 1;



				}
			}
			if (DSEn == 1)
			{
				l2count++;
				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
				next - > readFromAddressL2(addr);
		}

		if (invalid == 0) //when all are filled, check highest LRU_count and update tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = Arr[indexr][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L1_ASSOC; j++)
			{
				if (Arr[indexr][j].t_LRU_counter > temp)
				{
					temp = Arr[indexr][j].t_LRU_counter;
					highest_LRU = j;
				}
			}
			unsigned int z = (Arr[indexr][highest_LRU].t_tag << indexsize | indexr) << bo;

			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{

				Arr[indexr][k].t_LRU_counter++;
			}
			Arr[indexr][highest_LRU].t_tag = tagr;
			Arr[indexr][highest_LRU].t_LRU_counter = 0;
			if (Arr[indexr][highest_LRU].t_dirty_bit == 1)
			{
				write_backs++;
				if (DSEn == 1)
					next - > writeToAddressDS(z);
				if (L2En == 1)
					next - > writeToAddressL2(z);
			}
			Arr[indexr][highest_LRU].t_dirty_bit = 0;
			if (DSEn == 1)
			{
				l2count++;
				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
				next - > readFromAddressL2(addr);

		}

		if (invalid != 0 && invalid != 1) //case when multiple blocks are empty
		{
			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[indexr][k].t_valid_bit == 0 && Arr[indexr][k].t_LRU_counter > currentHighest) //check which are empty and find the
				{ //one with highest lru count
					currentHighest = Arr[indexr][k].t_LRU_counter;
					highest_LRU = k;
				}

			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[indexr][k].t_LRU_counter < Arr[indexr][highest_LRU].t_LRU_counter)
					Arr[indexr][k].t_LRU_counter++;
			}

			Arr[indexr][highest_LRU].t_tag = tagr; //update tag,lru counter and valid bit 
			Arr[indexr][highest_LRU].t_LRU_counter = 0;
			Arr[indexr][highest_LRU].t_valid_bit = 1;
			if (DSEn == 1)
			{
				l2count++;
				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
				next - > readFromAddressL2(addr);

		}

	}

}

void L1CACHE::writeToAddressL1(unsigned int addr)
{
	writes++;
	int w_hit_block = 0;
	unsigned int tag = addr >> (bo + indexsize); //extract tag from address
	unsigned int index = ((1 << indexsize) - 1) & (addr >> (bo)); //extract index from address
	unsigned int block_bits = ((1 << bo) - 1) & (addr >> 1);
	int w_hit = 0;
	for (unsigned int j = 0; j < L1_ASSOC; j++)
	{
		if (Arr[index][j].t_tag == tag && Arr[index][j].t_valid_bit == 1) //find if hit or miss and if hit find the hit block
		{
			w_hit = 1;
			w_hit_block = j;
		}
	}

	if (w_hit == 1) //hit case
	{
		unsigned int oldLRU;
		oldLRU = Arr[index][w_hit_block].t_LRU_counter;
		for (unsigned int i = 0; i < L1_ASSOC; i++)
		{
			if (Arr[index][i].t_LRU_counter < oldLRU) //update counters of other blocks which have lru counter value
				Arr[index][i].t_LRU_counter += 1; //less than hit block counter value
		}
		Arr[index][w_hit_block].t_LRU_counter = 0; //put lru counter value of hit block to zero and put dirty flag to 1
		Arr[index][w_hit_block].t_dirty_bit = 1;
	}


	if (w_hit == 0) //miss case
	{
		write_misses++;

		int invalid = 0;
		for (unsigned int vl = 0; vl < L1_ASSOC; vl++)
		{
			if (Arr[index][vl].t_valid_bit == 0) //check for number of invalids
			{
				invalid++;

			}
		}


		if (invalid == 1) //case when only one block is empty
		{
			if (DSEn == 1)
			{
				l2count++;

				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
				next - > readFromAddressL2(addr);
			for (unsigned int i = 0; i < L1_ASSOC; i++)
			{
				if (Arr[index][i].t_valid_bit == 0) //update tag,LRU_counter and valid_bit when only one block is invalid(empty)
				{
					for (unsigned int j = 0; j < L1_ASSOC; j++)
					{

						Arr[index][j].t_LRU_counter++;
					}

					Arr[index][i].t_tag = tag;
					Arr[index][i].t_LRU_counter = 0;
					Arr[index][i].t_valid_bit = 1;
					Arr[index][i].t_dirty_bit = 1;


				}
			}
		}

		if (invalid == 0) //when none are empty, check highest LRU_count and add tag to that block
		{
			unsigned int temp, highest_LRU;
			temp = Arr[index][0].t_LRU_counter;
			highest_LRU = 0;
			for (unsigned int j = 0; j < L1_ASSOC; j++)
			{
				if (Arr[index][j].t_LRU_counter > temp)
				{
					temp = Arr[index][j].t_LRU_counter;
					highest_LRU = j;
				}
			}
			unsigned int z = (Arr[index][highest_LRU].t_tag << indexsize | index) << bo;

			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{

				Arr[index][k].t_LRU_counter++;
			}

			Arr[index][highest_LRU].t_tag = tag;
			Arr[index][highest_LRU].t_LRU_counter = 0;
			if (Arr[index][highest_LRU].t_dirty_bit == 1)
			{
				write_backs++;
				if (DSEn == 1)
					next - > writeToAddressDS(z);
				if (L2En == 1)
					next - > writeToAddressL2(z);
			}
			if (DSEn == 1)
			{
				l2count++;
				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
				next - > readFromAddressL2(addr);
			Arr[index][highest_LRU].t_dirty_bit = 1;


		}

		if (invalid != 0 && invalid != 1) //when multiple are empty 
		{
			unsigned int currentHighest = 0, highest_LRU;
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[index][k].t_valid_bit == 0 && Arr[index][k].t_LRU_counter > currentHighest) //find the empty blocks and find the
				{ //one with highest lru counter value
					currentHighest = Arr[index][k].t_LRU_counter;
					highest_LRU = k;
				}

			}
			if (DSEn == 1)
			{
				//l2count++;
				next - > readFromAddressDS(addr);
			}
			if (L2En == 1)
			{
				next - > readFromAddressL2(addr);
			}

			//Update the LRU Counters
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[index][k].t_LRU_counter < Arr[index][highest_LRU].t_LRU_counter)
					Arr[index][k].t_LRU_counter++;
			}

			Arr[index][highest_LRU].t_tag = tag;
			Arr[index][highest_LRU].t_LRU_counter = 0;
			Arr[index][highest_LRU].t_valid_bit = 1;
			Arr[index][highest_LRU].t_dirty_bit = 1;
		}

	}

}



void L1CACHE::printStatsL1()
{
	cout << "===== L1 Contents =====" << endl;

	for (unsigned int i = 0; i < sets; i++)
	{
		cout << "set  " << dec << i << ": ";
		for (unsigned int j = 0; j < L1_ASSOC; j++)
		{
			for (unsigned int k = 0; k < L1_ASSOC; k++)
			{
				if (Arr[i][k].t_LRU_counter == j)
				{
					cout << " " << hex << Arr[i][k].t_tag;
					if (Arr[i][k].t_dirty_bit == 1)
					{
						cout << " D";
					}
					else
					{
						cout << " N";
					}
					cout << " ||";
				}
			}
		}
		cout << '\n';
	}


}

void L2CACHE::printStatsL2()
{
	if (L2En == 1)
	{
		cout << "===== L2 Contents =====" << endl;

		for (unsigned int i = 0; i < sets; i++)
		{
			cout << "set  " << dec << i << ": ";
			for (unsigned int j = 0; j < L2_ASSOC; j++)
			{
				for (unsigned int k = 0; k < L2_ASSOC; k++)
				{
					if (ArrL2[i][k].t_LRU_counter == j)
					{
						cout << " " << hex << ArrL2[i][k].t_tag;
						if (ArrL2[i][k].t_dirty_bit == 1)
						{
							cout << " D";
						}
						else
						{
							cout << " N";
						}
						cout << " ||";
					}
				}
			}
			cout << '\n';
		}
	}
	if (DSEn == 1)
	{
		cout << "===== L2 Address Array Contents =====" << endl;

		for (unsigned int i = 0; i < sets; i++)
		{
			cout << "set  " << dec << i << ": ";
			for (unsigned int k = 0; k < L2_ADDR_TAGS; k++)
			{
				cout << " " << hex << L2_tagstore[i].ds_tag[k];

			}
			cout << " ||";
			cout << '\n';
		}
		cout << '\n';
		cout << "===== L2 Data Array contents =====" << endl;
		for (unsigned int i = 0; i < sets; i++)
		{
			cout << "set  " << dec << i << ": ";
			for (unsigned int j = 0; j < L2_DATA_BLOCKS; j++)
			{
				cout << L2_datastore[i].ds_sel_bit[j] << ",";
				if (L2_datastore[i].ds_valid_bit[j] == 1)
					cout << "V,";
				else
					cout << "I,";


				if (L2_datastore[i].ds_dirty_bit[j] == 1)
				{
					cout << "D";
				}
				else
				{
					cout << "N";
				}

				cout << "   ";

			}
			cout << " ||";
			cout << '\n';
		}
	}
}





int main(int argc, char * argv[])
{
	char * trace_file;
	unsigned int c_BLOCKSIZE = atoi(argv[1]); //command line arguments
	unsigned int c_L1_SIZE = atoi(argv[2]);
	unsigned int c_L1_ASSOC = atoi(argv[3]);
	unsigned int c_L2_SIZE = atoi(argv[4]);
	unsigned int c_L2_ASSOC = atoi(argv[5]);
	unsigned int c_L2_DATA_BLOCKS = atoi(argv[6]);
	unsigned int c_L2_ADDR_TAGS = atoi(argv[7]);
	unsigned int c_DSEn = 0;
	unsigned int c_L2En = 0;
	if (c_L2_ADDR_TAGS > 1)
	{
		c_DSEn = 1;
	}
	if ((c_L2_ADDR_TAGS == 1) && (c_L2_SIZE > 0))
	{
		c_L2En = 1;
	}
	trace_file = argv[8];
	L2CACHE L2(c_BLOCKSIZE, c_L2_SIZE, c_L2_ASSOC, c_L2_DATA_BLOCKS, c_L2_ADDR_TAGS);
	L2CACHE * c_next = & L2;
	L1CACHE L1(c_BLOCKSIZE, c_L1_SIZE, c_L1_ASSOC, c_next, c_DSEn, c_L2En);
	ifstream infile(trace_file);
	char instr;
	unsigned int addr;
	while (infile >> instr >> hex >> addr) //read the instructions and address from trace file
	{
		if (instr == 'r')
			L1.readFromAddressL1(addr);
		else if (instr == 'w')
			L1.writeToAddressL1(addr);

	}
	double miss_rate = 0.0 f;
	unsigned int total_mem_traffic = 0;

	if (L1.reads + L1.writes == 0)
		miss_rate = 1.0 f;
	else
		miss_rate = 1.0 f * (L1.read_misses + L1.write_misses) /
		(L1.reads + L1.writes);

	total_mem_traffic = L1.read_misses + L1.write_misses + L1.write_backs;


	double miss_rateL2 = 0.0 f;
	unsigned int total_mem_trafficL2 = 0;

	if (L2.readsL2 + L2.writesL2 == 0)
		miss_rateL2 = 1.0 f;
	else
		miss_rateL2 = 1.0 f * L2.read_missesL2 / L2.readsL2;

	total_mem_trafficL2 = L2.read_missesL2 + L2.write_missesL2 + L2.write_backsL2;

	printf("===== Simulator configuration =====\n");
	printf("BLOCKSIZE:               %u\n", c_BLOCKSIZE);
	printf("L1_SIZE:                 %u\n", c_L1_SIZE);
	printf("L1_ASSOC:                %u\n", c_L1_ASSOC);
	printf("L2_SIZE:                 %u\n", c_L2_SIZE);
	printf("L2_ASSOC:                %u\n", c_L2_ASSOC);
	printf("L2_DATA_BLOCKS:          %u\n", c_L2_DATA_BLOCKS);
	printf("L2_ADDRESS_TAGS:         %u\n", c_L2_ADDR_TAGS);
	printf("trace_file:              %s\n", trace_file);
	cout << '\n';
	L1.printStatsL1();
	cout << '\n';
	if (c_L2_SIZE > 0)
	{
		L2.printStatsL2();
	}

	printf("\n===== Simulation Results =====\n");
	printf("a. number of L1 reads: %u\n", L1.reads);
	printf("b. number of L1 read misses: %u\n", L1.read_misses);
	printf("c. number of L1 writes: %u\n", L1.writes);
	printf("d. number of L1 write misses: %u\n", L1.write_misses);
	printf("e. L1 miss rate: %.4f\n", miss_rate);
	printf("f. number of writebacks from L1 memory: %u\n", L1.write_backs);
	if (c_L2_SIZE > 0)
	{
		printf("g. number of L2 reads: %u\n", L2.readsL2);
		printf("h. number of L2 read misses: %u\n", L2.read_missesL2);
		printf("i. number of L2 writes: %u\n", L2.writesL2);
		printf("j. number of L2 write misses: %u\n", L2.write_missesL2);
		if (c_L2_ADDR_TAGS > 1 || c_L2_DATA_BLOCKS > 1)
		{
			printf("k. number of L2 sector misses: %u\n", L2.sector_miss);
			printf("l. number of L2 cache block misses: %u\n", L2.block_miss);
			printf("m. L2 miss rate: %.4f\n", miss_rateL2);
			printf("n. number of writebacks from L2 memory: %u\n", L2.write_backsL2);
			printf("o. total memory traffic: %u\n", total_mem_trafficL2);
		}
		else if (c_L2_ADDR_TAGS == 1 && c_L2_SIZE > 0 && c_L2_DATA_BLOCKS == 1)
		{
			printf("k. L2 miss rate: %.4f\n", miss_rateL2);
			printf("l. number of writebacks from L2 memory: %u\n", L2.write_backsL2);
			printf("m. total memory traffic: %u\n", total_mem_trafficL2);
		}
	}
	else if (c_L2_SIZE == 0)
		printf("g. total memory traffic: %u\n", total_mem_traffic);

}
