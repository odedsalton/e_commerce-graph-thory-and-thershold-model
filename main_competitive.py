import hw1_part2 as sw
import hw1_part2_basic as st
def main():
    links_data = sw.ModelData('links_dataset.csv')
    removed_list = sw.competitive_part(links_data,1000)
    sw.write_file_competition(removed_list)
    ########################### basic solution ################
    removed_list = st.removing_highest_weight(links_data,1000)
    sw.write_file_competition(removed_list)
    return
if __name__ == '__main__':
    main()