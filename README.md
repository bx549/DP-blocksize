# DP-blocksize
A dynamic block size model for transaction fee-only cryptocurrencies

After all Bitcoins have been minted, user transaction fees will be the only source of revenue to support the network operations, that is to say, the creation of new blocks. A certain amount of network congestion/delay will be necessary in order to keep transaction fees high enough to attract what is currently called "mining" operations. However, too much congestion/delay and unreasonable fees will cause users to seek other payment options. We propose a dynamic block size model for the post-coinbase-reward scenario. The block size is adjusted so that revenue for "mining" operations is maximized while congestion and user delay is kept to a reasonable level.

This code implements the value iteration algorithm of Dynamic Programming to obatin an optimal policy for a dynamic block size.

