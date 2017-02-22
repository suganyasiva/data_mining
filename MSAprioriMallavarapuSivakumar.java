package dm.project1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.util.Combinations;

public class MSAprioriMallavarapuSivakumar {
	LinkedHashMap<Integer, Double> minSup = new LinkedHashMap<Integer, Double>();
	LinkedList<Integer> must_have = new LinkedList<Integer>();
	LinkedList<LinkedList<Integer>> notogether = new LinkedList<LinkedList<Integer>>();
	double sdc;
	ArrayList<ArrayList<Integer>> transactions = new ArrayList<ArrayList<Integer>>();
	LinkedHashMap<Integer, Double> supportCount = new LinkedHashMap<Integer, Double>();
	LinkedHashMap<LinkedList<Integer>, Double> allSupportCounts = new LinkedHashMap<LinkedList<Integer>, Double>();
	LinkedList<LinkedList<LinkedList<Integer>>> F = new LinkedList<LinkedList<LinkedList<Integer>>>();
	LinkedHashMap<LinkedList<Integer>, Double> tailCount = new LinkedHashMap<LinkedList<Integer>, Double>();
	
	public static void main(String[] args) {
		MSAprioriMallavarapuSivakumar m = new MSAprioriMallavarapuSivakumar();
		// read the trasactions, MIS values , sdc values
		if(args.length != 0){
			m.setMinSup(args[0]);
			m.setTransactions(args[1]);
		} else {
			System.out.println("Provide input and parameter files as command-line arguments");
			System.exit(0);
		}
		m.setSupportCount();
		//set support counts for the candidates
		m.setAllSupportCounts();
		//the main MS-Apriori function
		m.compute();
		m.outputData();
	}

	public void outputData() {
		/*
		 * Writes the output to a file: output-file.txt
		 * Output has : Frequent itemsets, counts and tailcounts
		 */
		try {
			PrintWriter writer = new PrintWriter(".\\data\\output-file.txt", "UTF-8");
			for (int i = 0; i < F.size(); i++) {
				LinkedList<LinkedList<Integer>> nItemset = F.get(i);
				int itemsetNumber = i + 1;
				if(!nItemset.isEmpty()){
					System.out.println("Frequent " + itemsetNumber + "-itemsets");
					writer.println("Frequent " + itemsetNumber + "-itemsets\n");
					for (int j = 0; j < nItemset.size(); j++) {
						LinkedList<Integer> frequentItem = nItemset.get(j);
						System.out
								.println("\t" + this.allSupportCounts.get(frequentItem).intValue() + " : " + frequentItem.toString().replace("[", "{").replace("]", "}"));
						writer.println("\t" + this.allSupportCounts.get(frequentItem).intValue() + " : " + frequentItem.toString().replace("[", "{").replace("]", "}"));
						if (i != 0) {
							System.out.println("Tailcount: " + this.tailCount.get(frequentItem).intValue());
							writer.println("Tailcount: " + this.tailCount.get(frequentItem).intValue());
						}
					}
					System.out.println("\n\tTotal number of Frequent " + itemsetNumber + "-itemsets = " + nItemset.size());
					writer.println("\n\tTotal number of Frequent " + itemsetNumber + "-itemsets = " + nItemset.size() + "\n\n");
				}
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void setAllSupportCounts() {
		/*
		 * initialises the hashmap of support counts with the support count of each item in Transaction
		 */
		Iterator<Entry<Integer, Double>> it = this.getSupportCount().entrySet().iterator();
		int tSize = this.transactions.size();
		while (it.hasNext()) {
			Entry<Integer, Double> entry = it.next();
			Integer item = entry.getKey();
			Double count = entry.getValue();
			this.allSupportCounts.put(new LinkedList<Integer>(Arrays.asList(item)), count * tSize);
		}
	}

	public LinkedHashMap<Integer, Double> sortHashMapByValues(HashMap<Integer, Double> minsuphash) {
		List<Double> mapValues = new ArrayList<Double>(minsuphash.values());
		List<Integer> mapKeys = new ArrayList<Integer>(minsuphash.keySet());
		Collections.sort(mapValues);
		Collections.sort(mapKeys);

		LinkedHashMap<Integer, Double> sortedMap = new LinkedHashMap<Integer, Double>();

		Iterator<Double> valueIt = mapValues.iterator();
		while (valueIt.hasNext()) {
			Double val = valueIt.next();
			Iterator<Integer> keyIt = mapKeys.iterator();

			while (keyIt.hasNext()) {
				Integer key = keyIt.next();
				Double comp1 = minsuphash.get(key);
				Double comp2 = val;

				if (comp1.equals(comp2)) {
					keyIt.remove();
					sortedMap.put(key, val);
					break;
				}
			}
		}

//		System.out.println("M: " + sortedMap.toString());
		this.minSup = sortedMap;
		return this.minSup;
	}

	public LinkedHashMap<Integer, Double> calc_support_count() {
		/*
		 * Calculates support count of 1 item set only
		 */
		LinkedHashMap<Integer, Double> support_counts = new LinkedHashMap<Integer, Double>();
		/* 
		 * populate it with the linkedhashmap of itemsets unique elements as
		 * in MIS and initialise to 0
		 */
		Iterator<Entry<Integer, Double>> it = this.getMinSup().entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			Integer key = (Integer) pair.getKey();
			support_counts.put(key, Double.parseDouble("0.0"));
		}

		ArrayList<ArrayList<Integer>> transactions = this.getTransac();
		Iterator<ArrayList<Integer>> outer = transactions.iterator();
		while (outer.hasNext()) {
			ArrayList<Integer> sub_transactions = outer.next();
			Iterator<Integer> inner = sub_transactions.iterator();
			while (inner.hasNext()) {
				Integer key = inner.next();
				double value = support_counts.get(key);
				support_counts.put(key, value + 1);
			}
		}
		it = support_counts.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			Integer key = (Integer) pair.getKey();
			Double value = (Double) pair.getValue();
			support_counts.put(key, value / this.transactions.size());
		}
		return support_counts;
	}

	public LinkedList<LinkedList<Integer>> genLevel2(LinkedList<Integer> L) {
		/*
		 * generates C2
		 */
		LinkedList<LinkedList<Integer>> c2 = new LinkedList<LinkedList<Integer>>();
		double local_sdc = this.getSdc();
		LinkedHashMap<Integer, Double> local_sup = this.getSupportCount();
		LinkedHashMap<Integer, Double> local_MIS = this.getMinSup();
		Iterator<Integer> it_L = L.iterator();
		for (int i = 0; i < L.size(); i++) {
			Integer key = L.get(i);
			if (local_sup.get(key) >= local_MIS.get(key)) {
				for (int j = i + 1; j < L.size(); j++) {
					Integer sub_key = L.get(j);
					double diff = Math.abs(local_sup.get(sub_key) - local_sup.get(key));
					if ((local_sup.get(sub_key) >= local_MIS.get(key)) && diff <= local_sdc) {
						LinkedList<Integer> tempc = new LinkedList<Integer>();
						tempc.add(key);
						tempc.add(sub_key);
						c2.add(tempc);
					}
				}
			}
		}
//		System.out.println("C2: " + c2);
		return c2;
	}

	public boolean findinnothave(LinkedList<Integer> fkToCheck) {
		/*
		 * checks if the elements in not together are
		 * there in fk or not returns false if yes must take input as a
		 * linkedlist of integers
		  */
		LinkedList<LinkedList<Integer>> local_notogether = new LinkedList<LinkedList<Integer>>();
		local_notogether = this.getNotogether(); 
		int count = 0;
		for (LinkedList<Integer> mustNotHave : local_notogether) {
			if (fkToCheck.containsAll(mustNotHave)) {
				return false;
			}
		}
		return true;
	}

	public boolean containsAny(LinkedList<Integer> itemSet) {
		LinkedList<Integer> mustHave = this.getMust_have();
		for (int i = 0; i < mustHave.size(); i++) {
			if (itemSet.contains(mustHave.get(i))) {
				return true;
			}
		}
		return false;
	}

	public void compute() {
		/*
		 * The MS-Apriori algorithm
		 */
		LinkedList<Integer> L = this.init_pass();
		LinkedList<LinkedList<Integer>> F1 = this.getF1(L);
		LinkedList<LinkedList<LinkedList<Integer>>> C = new LinkedList<LinkedList<LinkedList<Integer>>>();
		LinkedHashMap<LinkedList<Integer>, Double> candidateCount = new LinkedHashMap<LinkedList<Integer>, Double>();
		F.add(F1);
		//current is used to know F(k-1) for Ck
		LinkedList<LinkedList<Integer>> current = new LinkedList<LinkedList<Integer>>();
		current = F1;
		int k = 2; // count the number of F generated
		while (k >= 2) {
			if (current.isEmpty()) {
				break;
			} else {
				// present is used to know ck for fk
				LinkedList<LinkedList<Integer>> present = new LinkedList<LinkedList<Integer>>();
				if (k == 2) {
					/*
					 * Level2 generation function
					 */
					LinkedList<LinkedList<Integer>> c2 = this.genLevel2(L);
					present = c2;
					C.add(c2);
				} else {
					LinkedList<LinkedList<Integer>> ck = this.MSCandidateGen(current);
					present = ck;
					C.add(ck);
				}
				// for computing support counts of candidates, initialise to 0.0
				for (LinkedList<Integer> candidate : present) {
					candidateCount.put(candidate, Double.parseDouble("0.0"));
				}
				for (LinkedList<Integer> candidate : present) {
					LinkedList<Integer> tail = new LinkedList<Integer>();
					tail.addAll(candidate);
					tail.removeFirst(); // removes first element initailises the tail count hashmap
					this.tailCount.put(candidate, Double.parseDouble("0.0"));
				}
				ArrayList<ArrayList<Integer>> local_transaction = this.getTransac();
				// look for candidate in each transaction
				for (ArrayList<Integer> transaction : local_transaction) {
					for (LinkedList<Integer> candidate : present) {
						if (containsCandidate(transaction, candidate)) {
							candidateCount.replace(candidate, candidateCount.get(candidate) + 1);
						}
						// add new tailcounts for the candidates
						LinkedList<Integer> tail = new LinkedList<Integer>();
						tail.addAll(candidate);
						tail.removeFirst();
						if (containsCandidate(transaction, tail)) {
							this.tailCount.replace(candidate, this.tailCount.get(candidate) + 1); 
						}
					}
				}

				k++;
				this.allSupportCounts.putAll(candidateCount);
				LinkedList<LinkedList<Integer>> Fk = this.getFk(present, candidateCount);
				current = Fk; // marks f(k-1) for ck
				if (Fk.size() > 0) {
					F.add(Fk);
				}
			}
		}

		/* Check for must-have and cannot-be-together */
		LinkedList<LinkedList<LinkedList<Integer>>> FCopy = new LinkedList<LinkedList<LinkedList<Integer>>>();
		for (LinkedList<LinkedList<Integer>> tempsub : F) {
			FCopy.add(new LinkedList<LinkedList<Integer>>(tempsub));
		}

		for (int i = 0; i < FCopy.size(); i++) {
			for (int j = 0; j < FCopy.get(i).size(); j++) {
				if (!containsAny(FCopy.get(i).get(j)) || !findinnothave(FCopy.get(i).get(j))) {
					F.get(i).remove(FCopy.get(i).get(j));
				}
			}
		}

//		System.out.println("F: " + this.F);
	}

	public boolean containsCandidate(ArrayList<Integer> transaction, LinkedList<Integer> candidate) {
		/* function to check if a candidate is part of a given transaction */
		int count = 0;
		for (Integer item : candidate) {
			if (transaction.contains(item)) {
				count++;
			}
		}
		return (count == candidate.size()) ? true : false;
	}

	public LinkedList<LinkedList<Integer>> MSCandidateGen(LinkedList<LinkedList<Integer>> Fk) {
		/* function to generate candidates for level > 2 */
		LinkedList<LinkedList<Integer>> ck = new LinkedList<LinkedList<Integer>>();
		LinkedHashMap<Integer, Double> local_MIS = this.getMinSup();
		LinkedList<LinkedList<Integer>> Fksorted = new LinkedList<LinkedList<Integer>>(Fk);
		for (int j = 0; j < Fksorted.size(); j++) {
			LinkedList<Integer> itemset = new LinkedList<Integer>();
			itemset = Fksorted.get(j);
			for (int i = j + 1; i < Fksorted.size(); i++) {
				LinkedList<Integer> c = new LinkedList<Integer>();
				LinkedList<LinkedList<Integer>> subsetsOfC = new LinkedList<LinkedList<Integer>>();
				LinkedList<Integer> f2 = Fksorted.get(i);
				Integer lastItemInF1 = itemset.getLast();
				Integer lastItemInF2 = f2.getLast();
				if (differsOnlyInLastItem(itemset, f2) && isLessThan(local_MIS.get(lastItemInF1), local_MIS.get(lastItemInF2)) && (Math
						.abs(this.supportCount.get(lastItemInF1) - this.supportCount.get(lastItemInF2)) <= this.sdc)) {
						c.addAll(itemset);
						c.add(lastItemInF2);
						ck.add(c);
					subsetsOfC = findKMinusOneSubsetsOfC(c);
					for (LinkedList<Integer> subset : subsetsOfC) {
						if (subset.contains(c.getFirst())
								|| this.minSup.get(c.getFirst()) == this.minSup.get(c.get(1))) {
							if (!Fksorted.contains(subset)) { 
								ck.remove(c);
							}
						}
					}
				}
			}
		}
		LinkedList<LinkedList<Integer>> unsortck = new LinkedList<LinkedList<Integer>>();
		unsortck.addAll(ck);
		return unsortck;
	}

	public LinkedList<LinkedList<Integer>> findKMinusOneSubsetsOfC(LinkedList<Integer> c) {
		/*
		 * finds the subsets of ck
		 */
		LinkedList<LinkedList<Integer>> subSetsOfC = new LinkedList<LinkedList<Integer>>();
		for (Iterator<int[]> iter = new Combinations(c.size(), c.size() - 1).iterator(); iter.hasNext();) {
			LinkedList<Integer> subset = new LinkedList<Integer>();
			for (int i : iter.next()) {
				subset.add(c.get(i));
			}
			subSetsOfC.add(subset);
		}
		return subSetsOfC;
	}

	public boolean differsOnlyInLastItem(LinkedList<Integer> f1, LinkedList<Integer> f2) {
		/*
		 * Find the f1 f2 to be joined for Fk
		 */
		LinkedList<Integer> list1 = new LinkedList<Integer>();
		LinkedList<Integer> list2 = new LinkedList<Integer>();
		list1.addAll(f1);
		list2.addAll(f2);
		if (!list1.removeLast().equals(list2.removeLast()) && list1.equals(list2)) {
			return true;
		}
		return false;
	}

	public boolean isLessThan(Double integer1, Double integer2) {
		return (integer1 >= integer2) ? false : true;
	}

	public LinkedList<LinkedList<Integer>> getFk(LinkedList<LinkedList<Integer>> ck,
			LinkedHashMap<LinkedList<Integer>, Double> candidateCount) {
		/*
		 * function to generate frequent itemsets candidateCount -> counts of
		 * each subset (still a set) of ck
		 */
		LinkedList<LinkedList<Integer>> Fk = new LinkedList<LinkedList<Integer>>();
		LinkedHashMap<Integer, Double> local_MIS = this.getMinSup();

		for (LinkedList<Integer> c : ck) {
			Iterator<Integer> it = c.iterator();
			double firstMIS = local_MIS.get(it.next()); // gets MIS of c[1]
			if (candidateCount.get(c) / this.transactions.size() >= firstMIS) {
				Fk.add(c);
			}
		}
//		System.out.println("Fk: " + Fk);
		return Fk;
	}

	public LinkedList<LinkedList<Integer>> getF1(LinkedList<Integer> L) {
		/* function to generate level 1 frequent itemset */

		// array List to be returned
		LinkedList<LinkedList<Integer>> F1 = new LinkedList<LinkedList<Integer>>();
		// support counts (computed already do not divide by n again)
		LinkedHashMap<Integer, Double> sup_count = this.getSupportCount();
		// MIS values for each element in itemset
		LinkedHashMap<Integer, Double> local_minsup = this.getMinSup();
		Iterator<Integer> it = L.iterator();
		// compare the MIS(i) to sup(i)
		while (it.hasNext()) {
			Integer key = it.next();
			double mis_key = local_minsup.get(key);
			double sup = sup_count.get(key);
			if (sup >= mis_key) {
				LinkedList<Integer> sub_set = new LinkedList<Integer>();
				sub_set.add(key);
				F1.add(sub_set);
			}
		}
//		System.out.println("F1: " + F1);
		return F1;
	}

	public LinkedList<Integer> init_pass() {
		/*
		 * Init pass step: 
		 * 1) Sorts the MINSUP hash according to the MINSUP
		 * 2) return L (candidates for 1-set frequent items
		 */
		// sorts the MIS for lowest first
		LinkedHashMap<Integer, Double> M = this.sortHashMapByValues(this.minSup);
		// gets the support counts for 1- itemset
		LinkedHashMap<Integer, Double> support_counts = this.calc_support_count();
		// these are all the transactions
		ArrayList<ArrayList<Integer>> transac = this.getTransac();
		// List to be returned
		LinkedList<Integer> L = new LinkedList<Integer>();
		Iterator<Entry<Integer, Double>> it_support_count = support_counts.entrySet().iterator();
		Iterator<Entry<Integer, Double>> it = M.entrySet().iterator();
		double mis_min = 0.0;
		int position = 0;
		while (it.hasNext()) {
			Map.Entry pair_sorted = (Map.Entry) it.next();
			mis_min = (Double) pair_sorted.getValue();
			double sup_min = (Double) support_counts.get(pair_sorted.getKey());
			if (sup_min >= mis_min) {
				L.add((Integer) pair_sorted.getKey());
				break;
			}
			position++;
		}
		int pointer = 0;
		while (it_support_count.hasNext() && pointer <= position) {
			// compare for subsequent elements only not all
			it_support_count.next();
			pointer++;
		}
		while (it_support_count.hasNext()) {
			Map.Entry pair_sup = (Map.Entry) it_support_count.next();
			Integer key_sup = (Integer) pair_sup.getKey();
			double sup = (Double) pair_sup.getValue();
			if (sup >= mis_min) {
				// support is greater than minimum support include it
				L.add(key_sup);
			}
		}
//		System.out.println("L:" + L);
		return L;
	}

	public void setSdc(String sdc) {
		/*
		 * reads the SDC into the variable
		 */
		try {
			this.sdc = Double.parseDouble(sdc.split("=")[1]);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void setMinSup(String filename) {
		/*
		 * Read minimum support constraints from file return HashMap of MIS
		 * values needs to be sorted
		 */
		try {
			String path= ".\\data\\".concat(filename);
			BufferedReader br = new BufferedReader(new FileReader(path));
			String next = new String();
			while ((next = br.readLine()) != null) {
				if (next.contains("MIS")) {
					String split[] = next.split("=");
					Matcher m = Pattern.compile("\\(([^)]+)\\)").matcher(split[0]);
					while (m.find()) {
						this.minSup.put(Integer.parseInt(m.group(1)),
								Double.parseDouble(split[split.length - 1].trim()));
					}
				} else if (next.contains("SDC")) {
					this.setSdc(next);
				} else if (next.contains("cannot_be_together:")) {
					this.setnotogether(next);
				} else if (next.contains("must-have:")) {
					this.setmustHave(next);
				}
			}
			br.close();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public void setTransactions(String filename) {
		/*
		 * reads the transactions from the file into the variable
		 */
		try {
			String path= ".\\data\\".concat(filename);
			BufferedReader br = new BufferedReader(new FileReader(path));
			String next = new String();
			while ((next = br.readLine()) != null) {
				ArrayList eachTransac = new ArrayList<Integer>();
				String split[] = next.split(","); // separate elements
				String temp = split[0].split("\\{")[1];
				split[0] = temp;
				temp = split[split.length - 1].split("\\}")[0];
				split[split.length - 1] = temp;
				for (int i = 0; i < split.length; i++) {
					eachTransac.add(Integer.parseInt(split[i].trim()));
				}
				this.transactions.add(eachTransac);
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public LinkedHashMap<Integer, Double> getSupportCount() {
		return this.supportCount;
	}

	public void setSupportCount() {
		this.supportCount = this.calc_support_count();
	}

	public LinkedHashMap<Integer, Double> getMinSup() {
		return this.minSup;
	}

	public double getSdc() {
		return this.sdc;
	}

	public ArrayList<ArrayList<Integer>> getTransac() {
		return this.transactions;
	}

	public void setmustHave(String next) {
		/*
		 * reads the must haves into the variable
		 */
		LinkedList<Integer> keep = new LinkedList<Integer>();
		String split = next.split("must-have:")[1];
		String parts[] = split.split(" or ");
		for (int i = 0; i < parts.length; i++) {
			keep.add(Integer.parseInt(parts[i].trim()));
//			System.out.println("must_have:" + parts[i]);
		}
		this.must_have = keep;
	}

	public LinkedList<LinkedList<Integer>> getNotogether() {
		return notogether;
	}

	public LinkedList<Integer> getMust_have() {
		return must_have;
	}

	public LinkedHashMap<LinkedList<Integer>, Double> getAllSupportCounts() {
		return allSupportCounts;
	}

	public LinkedHashMap<LinkedList<Integer>, Double> getTailCount() {
		return tailCount;
	}

	public void setTailCount(LinkedHashMap<LinkedList<Integer>, Double> tailCount) {
		this.tailCount = tailCount;
	}

	public void setnotogether(String next) {
		/*
		 * read the not together into the variable
		 */
		LinkedList<LinkedList<Integer>> eliminate = new LinkedList<LinkedList<Integer>>();
		Matcher m = Pattern.compile("\\{(.*?)\\}").matcher(next);
		LinkedList<String> temporary = new LinkedList<String>();
		while (m.find()) {
			for (int count = 1; count <= m.groupCount(); count++) {
				temporary.add(m.group(count));
			}
		}
		for (int i = 0; i < temporary.size(); i++) {
			String split[] = temporary.get(i).split(",");
			LinkedList<Integer> subeliminate = new LinkedList<Integer>();
			for (int count = 0; count < split.length; count++) {
				subeliminate.add(Integer.parseInt(split[count].trim()));
//				System.out.println("dont have:" + split[count]);
			}
			eliminate.add(subeliminate);
		}
		// set the global variable
		this.notogether = eliminate;
	}

	public int randomWithRange(int min, int max) {
		int range = (max - min) + 1;
		return (int) (Math.random() * range) + min;
	}

	public double randomWithRangedouble(double min, double max) {
		double range = (max - min);
		return (Math.random() * range) + min;
	}
}
